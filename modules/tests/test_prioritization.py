import shutil

import pytest

from modules.prioritization import (
    count_primary_reads, count_variants, prioritize_variants, private_variants,
    read_coverage, read_merged_diff, summarize_samples, trimmomatic_surviving_reads,
    variant_key, vcf_allele_key, write_mt_read_fraction
)

VCF_HEADER = (
    "##fileformat=VCFv4.0\n"
    '##INFO=<ID=AC,Number=1,Type=Integer,Description="">\n'
    '##INFO=<ID=AN,Number=1,Type=Integer,Description="">\n'
    '##INFO=<ID=Locus,Number=1,Type=String,Description="">\n'
    '##INFO=<ID=NtVarH,Number=1,Type=Float,Description="">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n"
)


def write(path, text):
    path.write_text(text)
    return str(path)


def test_variant_key():
    assert variant_key("3594C") == ("snp", 3594, "C")
    assert variant_key("310C(Y)") == ("snp", 310, "C")
    assert variant_key("310.C") == ("ins", 310, "C")
    assert variant_key("3107d") == ("del", 3107, 3107)
    assert variant_key("523-524d") == ("del", 523, 524)
    assert variant_key("weird") is None


def test_vcf_allele_key_matches_classifier_notation():
    assert vcf_allele_key(3594, "C", "T") == variant_key("3594T")
    # insertion of C after position 310
    assert vcf_allele_key(310, "T", "TC") == variant_key("310.C")
    # deletion of positions 523-524
    assert vcf_allele_key(522, "CAC", "C") == variant_key("523-524d")
    assert vcf_allele_key(714, "AC", "A") == variant_key("715d")


def test_private_variants(tmp_path):
    merged_diff = write(tmp_path / "s_merged_diff.csv",
                        ",RSRS,MHCS,rCRS\n146T,yes,,\n310C(Y),yes,yes,yes\n"
                        "709A,yes,,yes\n8935T,yes,yes,yes\n152C(Y),,yes,yes\n")
    rows = read_merged_diff(merged_diff)
    assert rows[0] == ("146T", True, False, False)
    assert private_variants(rows) == ["310C(Y)", "8935T"]


def test_prioritize_variants(tmp_path):
    ann = write(tmp_path / "ann.vcf", VCF_HEADER.format(samples="s1\ts2") +
                "chrM\t310\t.\tT\tC\t.\tPASS\tAC=1;AN=2;Locus=MT-DLOOP;NtVarH=0.5\tGT\t1\t0\n"
                "chrM\t8935\t.\tC\tT\t.\tPASS\tAC=2;AN=2;Locus=MT-ATP8;NtVarH=0.01\tGT\t1\t1\n")
    d1 = write(tmp_path / "s1.csv", ",RSRS,MHCS,rCRS\n310C(Y),yes,yes,yes\n8935T,yes,yes,yes\n")
    d2 = write(tmp_path / "s2.csv", ",RSRS,MHCS,rCRS\n8935T,yes,yes,yes\n9999A,yes,yes,yes\n146T,yes,,\n")
    out = tmp_path / "prio.txt"
    counts = prioritize_variants([("s1", d1), ("s2", d2)], ann, str(out))
    assert counts == {"s1": 2, "s2": 2}
    lines = [l.split("\t") for l in out.read_text().splitlines()]
    assert lines[0] == ["Variant Allele", "Sample", "Locus", "NtVarH"]
    # sorted by ascending NtVarH, unannotated variants last
    assert lines[1] == ["8935T", "s1,s2", "MT-ATP8", "0.01"]
    assert lines[2] == ["310C(Y)", "s1", "MT-DLOOP", "0.5"]
    assert lines[3] == ["9999A", "s2", "", ""]


def test_count_variants_and_summary(tmp_path):
    vcf = write(tmp_path / "s1.vcf", VCF_HEADER.format(samples="s1") +
                "chrM\t73\t.\tA\tG\t.\tPASS\tAC=1;AN=1\tGT:DP:HF\t1:100:1.0\n"
                "chrM\t150\t.\tC\tT,A\t.\tPASS\tAC=1,1;AN=1\tGT:DP:HF\t1:100:0.85,0.1\n"
                "chrM\t200\t.\tA\tG\t.\tPASS\tAC=1;AN=1\tGT:DP:HF\t1:100:0.3\n"
                "chrM\t300\t.\tT\tC\t.\tPASS\tAC=1;AN=1\tGT:DP:HF\t1:1000:0.985\n"
                "chrM\t400\t.\tG\tA\t.\tPASS\tAC=1;AN=1\tGT:DP:HF\t1:1000:0.01\n")
    # homoplasmic: 1.0, 0.985; heteroplasmic: 0.85, 0.1, 0.3; low-level: 0.01
    assert count_variants(vcf, homoplasmy_threshold=0.97, heteroplasmy_min=0.03) == (6, 2, 3, 1)

    cov = write(tmp_path / "s1.cov", "chrM\t1\t10\nchrM\t2\t2\nchrM\t3\t0\nchrM\t4\t20\n")
    assert read_coverage(cov, min_depth=5) == (50.0, 8.0)

    best = write(tmp_path / "s1_best.csv", "s1,H1;H1a\n")
    mt_reads = write(tmp_path / "s1_mt.tsv", "sample\treads_after_trimming\tmtDNA_reads\tmtDNA_reads_pct\n"
                     "s1\t2000\t50\t2.5\n")
    out = tmp_path / "summary.txt"
    summarize_samples([{"sample": "s1", "best_results": best, "vcf": vcf, "coverage": cov,
                        "mt_reads": mt_reads}],
                      {"s1": 3}, str(out), homoplasmy_threshold=0.97,
                      heteroplasmy_min=0.03, min_depth=5)
    row = out.read_text().splitlines()[-1].split("\t")
    assert row == ["s1", "2000", "50", "2.5", "50.0", "8.0", "H1;H1a", "6", "2", "3", "1", "3"]


TRIMMOMATIC_LOG = ("TrimmomaticPE: Started with arguments:\n ...\n"
                   "Input Read Pairs: 1000 Both Surviving: 900 (90.00%) Forward Only Surviving: 40 (4.00%) "
                   "Reverse Only Surviving: 10 (1.00%) Dropped: 50 (5.00%)\n"
                   "TrimmomaticPE: Completed successfully\n")

SAM = ("@SQ\tSN:chrM\tLN:16569\n"
       "r1\t99\tchrM\t100\t60\t10M\t=\t150\t60\tACGTACGTAC\tIIIIIIIIII\n"
       "r1\t147\tchrM\t150\t60\t10M\t=\t100\t-60\tACGTACGTAC\tIIIIIIIIII\n"
       "r2\t0\tchrM\t200\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
       "r2\t2048\tchrM\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
       "r3\t256\tchrM\t300\t0\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n")


def test_trimmomatic_surviving_reads(tmp_path):
    log = write(tmp_path / "trim.log", TRIMMOMATIC_LOG)
    # both mates of 900 pairs, plus 40 + 10 unpaired reads
    assert trimmomatic_surviving_reads(log) == 1850


@pytest.mark.skipif(shutil.which("samtools") is None, reason="needs samtools")
def test_mt_read_fraction(tmp_path):
    sam = write(tmp_path / "mt.sam", SAM)
    # primary reads only: not the supplementary (2048) and secondary (256) ones
    assert count_primary_reads(sam) == 3
    log = write(tmp_path / "trim.log", TRIMMOMATIC_LOG)
    out = tmp_path / "mt.tsv"
    write_mt_read_fraction("s1", [log, log], [sam, sam], str(out))
    assert out.read_text().splitlines()[1].split("\t") == ["s1", "3700", "6", "0.1622"]


def test_annotate_vcf_lock_free(tmp_path):
    # mtoolnote's database opened read-only without locks gives the same
    # annotation as mtoolnote's default connection
    mtoolnote = pytest.importorskip("mtoolnote")
    from modules.annotation import annotate_vcf
    vcf = write(tmp_path / "in.vcf", "##fileformat=VCFv4.0\n"
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="">\n'
                '##INFO=<ID=AC,Number=1,Type=Integer,Description="">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
                "chrM\t3243\t.\tA\tG\t.\tPASS\tAC=1\tGT\t1\n")
    annotate_vcf(vcf, str(tmp_path / "out.vcf"), "human")
    record = [l for l in (tmp_path / "out.vcf").read_text().splitlines() if not l.startswith("#")][0]
    assert "Locus=MT-TL1" in record
