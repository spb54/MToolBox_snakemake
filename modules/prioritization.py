#!/usr/bin/env python
""" Variant prioritization and per-sample summary for human mtDNA.

Port of MToolBox v1 prioritization.py/summary.py. A variant is
prioritized when it is private to the sample, i.e. it differs from all
three of RSRS, the MHCS of the predicted haplogroup and rCRS (the three
'yes' columns in the classifier _merged_diff.csv). Prioritized variants
are annotated with the mtoolnote INFO fields of the annotated VCF and
sorted by ascending nucleotide variability (NtVarH).
"""
import csv
import gzip
import math
import re
import subprocess
from collections import OrderedDict

SNP_RE = re.compile(r"^(\d+)([A-Z])(?:\(([A-Z])\))?$")
INS_RE = re.compile(r"^(\d+)\.([A-Z]+)$")
DEL_RE = re.compile(r"^(\d+)(?:-(\d+))?d$")

# INFO fields written by mtVariantCaller, not by mtoolnote
BASE_INFO = ("AC", "AN")
MISSING = {"", ".", "NA", "nan", "None"}


def variant_key(variant):
    """ Return a hashable key for a classifier variant notation.

    Examples: '3594C' and '310C(Y)' -> ('snp', 3594, 'C'),
    '310.C' -> ('ins', 310, 'C'), '523-524d' -> ('del', 523, 524).
    Returns None if the notation is not recognized.
    """
    m = SNP_RE.match(variant)
    if m:
        return ("snp", int(m.group(1)), m.group(2))
    m = INS_RE.match(variant)
    if m:
        return ("ins", int(m.group(1)), m.group(2))
    m = DEL_RE.match(variant)
    if m:
        start = int(m.group(1))
        end = int(m.group(2)) if m.group(2) else start
        return ("del", start, end)
    return None


def vcf_allele_key(pos, ref, alt):
    """ Return the variant_key() of a VCF allele (REF/ALT at POS). """
    if len(ref) == 1 and len(alt) == 1:
        return ("snp", pos, alt)
    if len(alt) > len(ref) and alt.startswith(ref):
        # insertion after the last reference base
        return ("ins", pos + len(ref) - 1, alt[len(ref):])
    if len(ref) > len(alt) and ref.startswith(alt):
        return ("del", pos + len(alt), pos + len(ref) - 1)
    return ("other", pos, ref, alt)


def _open(path):
    with open(path, "rb") as fh:
        gzipped = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gzipped else open(path)


def read_vcf(path):
    """ Parse a (possibly bgzipped) VCF.

    Returns:
        (info_ids, samples, records): INFO ids in header order, sample
        names, and one dict per ALT allele with keys pos, ref, alt,
        key, info (dict) and calls (dict sample -> dict of FORMAT
        fields, restricted to this allele where the field has one value
        per allele).
    """
    info_ids, samples, records = [], [], []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##INFO=<ID="):
                info_ids.append(line[len("##INFO=<ID="):].split(",")[0])
                continue
            if line.startswith("#CHROM"):
                samples = line.split("\t")[9:]
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            pos, ref, alts = int(fields[1]), fields[3], fields[4].split(",")
            info = OrderedDict()
            for item in fields[7].split(";"):
                if item:
                    k, _, v = item.partition("=")
                    info[k] = v
            fmt = fields[8].split(":") if len(fields) > 8 else []
            for i, alt in enumerate(alts):
                allele_info = OrderedDict(
                    (k, _split_allele(v, i, len(alts))) for k, v in info.items())
                calls = {}
                for name, value in zip(samples, fields[9:]):
                    call = dict(zip(fmt, value.split(":")))
                    calls[name] = {k: _split_allele(v, i, len(alts))
                                   for k, v in call.items()}
                records.append({"pos": pos, "ref": ref, "alt": alt,
                                "key": vcf_allele_key(pos, ref, alt),
                                "info": allele_info, "calls": calls})
    return info_ids, samples, records


def _split_allele(value, i, n_alts):
    """ Pick the i-th of n_alts comma separated values, if split per allele. """
    parts = value.split(",")
    return parts[i] if n_alts > 1 and len(parts) == n_alts else value


def read_merged_diff(path):
    """ Parse a classifier _merged_diff.csv.

    Returns:
        list of (variant, in_rsrs, in_mhcs, in_rcrs)
    """
    rows = []
    with open(path) as fh:
        reader = csv.reader(fh)
        next(reader, None)
        for row in reader:
            if row:
                row = row + [""] * (4 - len(row))
                rows.append((row[0],) + tuple(x == "yes" for x in row[1:4]))
    return rows


def private_variants(merged_diff_rows):
    """ Variants differing from RSRS, MHCS and rCRS alike. """
    return [r[0] for r in merged_diff_rows if r[1] and r[2] and r[3]]


def read_best_haplogroups(path):
    """ Return the best haplogroup(s) from a per-sample best results file. """
    with open(path) as fh:
        rows = [r for r in csv.reader(fh) if r]
    return rows[-1][1] if rows else ""


def read_coverage(path, min_depth=5):
    """ Return (percent of positions with depth >= min_depth, mean depth)
    from `samtools depth -a` output. """
    n = covered = total = 0
    with open(path) as fh:
        for line in fh:
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            depth = int(fields[2])
            n += 1
            total += depth
            covered += depth >= min_depth
    if n == 0:
        return 0.0, 0.0
    return round(100.0 * covered / n, 2), round(float(total) / n, 2)


def trimmomatic_surviving_reads(log_file):
    """ Number of reads left after trimming, from a Trimmomatic PE log:
    both mates of surviving pairs plus surviving unpaired reads. """
    with open(log_file) as fh:
        text = fh.read()
    m = re.search(r"Input Read Pairs: \d+ Both Surviving: (\d+) .*?"
                  r"Forward Only Surviving: (\d+) .*?Reverse Only Surviving: (\d+)", text)
    if m is None:
        raise ValueError("No Trimmomatic PE summary found in {}".format(log_file))
    both, forward, reverse = (int(x) for x in m.groups())
    return 2 * both + forward + reverse


def count_primary_reads(bam_file, samtools="samtools"):
    """ Number of primary mapped reads in a BAM (no secondary or
    supplementary alignments, e.g. of reads crossing the mtDNA origin). """
    out = subprocess.run([samtools, "view", "-c", "-F", "0x904", bam_file],
                         check=True, stdout=subprocess.PIPE,
                         universal_newlines=True).stdout
    return int(out.strip())


def write_mt_read_fraction(sample, trimmomatic_logs, mt_bams, out_file):
    """ Write the number and percentage of reads mapped to the mtDNA.

    Args:
        sample: sample name
        trimmomatic_logs: Trimmomatic logs of the sample's libraries
        mt_bams: mtDNA alignments of the sample's libraries, before
            duplicate removal (the total also includes duplicates)
        out_file: output tsv path
    """
    total = sum(trimmomatic_surviving_reads(f) for f in trimmomatic_logs)
    mt_reads = sum(count_primary_reads(f) for f in mt_bams)
    pct = round(100.0 * mt_reads / total, 4) if total else 0.0
    with open(out_file, "w") as out:
        out.write("sample\treads_after_trimming\tmtDNA_reads\tmtDNA_reads_pct\n")
        out.write("{}\t{}\t{}\t{}\n".format(sample, total, mt_reads, pct))


def read_mt_read_fraction(path):
    """ Return (reads_after_trimming, mtDNA_reads, mtDNA_reads_pct) as strings. """
    with open(path) as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r]
    return tuple(rows[-1][1:4])


def _to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _nt_variability(annotations):
    value = _to_float(annotations.get("NtVarH"))
    return math.inf if value is None else value


def prioritize_variants(sample_merged_diffs, annotated_vcf, out_file):
    """ Write prioritized variants of all samples, one row per variant.

    Args:
        sample_merged_diffs: list of (sample, merged_diff.csv path)
        annotated_vcf: mtoolnote annotated VCF with all samples
        out_file: output tsv path
    Returns:
        dict sample -> number of prioritized variants
    """
    info_ids, _, records = read_vcf(annotated_vcf)
    ann_ids = [i for i in info_ids if i not in BASE_INFO]
    annotations = {}
    for rec in records:
        annotations.setdefault(rec["key"], rec["info"])

    variants = OrderedDict()
    counts = {}
    for sample, merged_diff in sample_merged_diffs:
        private = private_variants(read_merged_diff(merged_diff))
        counts[sample] = len(private)
        for variant in private:
            variants.setdefault(variant, []).append(sample)

    rows = []
    for variant, samples in variants.items():
        ann = annotations.get(variant_key(variant), {})
        values = [ann.get(i, "") for i in ann_ids]
        values = ["" if v in MISSING else v for v in values]
        rows.append((_nt_variability(ann), variant,
                     ",".join(sorted(set(samples))), values))
    rows.sort(key=lambda r: (r[0], variant_key(r[1]) or ()))

    with open(out_file, "w") as out:
        out.write("\t".join(["Variant Allele", "Sample"] + ann_ids) + "\n")
        for _, variant, samples, values in rows:
            out.write("\t".join([variant, samples] + values) + "\n")
    return counts


def count_variants(sample_vcf, homoplasmy_threshold=0.97, heteroplasmy_min=0.03):
    """ Count the ALT alleles called in a single-sample VCF.

    Returns:
        (n_variants, n_homoplasmic, n_heteroplasmic, n_low_level):
        homoplasmic means HF >= homoplasmy_threshold, heteroplasmic
        heteroplasmy_min <= HF < homoplasmy_threshold and low-level
        HF < heteroplasmy_min
    """
    _, samples, records = read_vcf(sample_vcf)
    n = homo = het = low = 0
    for rec in records:
        call = rec["calls"].get(samples[0], {}) if samples else {}
        if call.get("GT", "1") in ("0", "."):
            continue
        n += 1
        hf = _to_float(call.get("HF"))
        if hf is None:
            continue
        if hf >= homoplasmy_threshold:
            homo += 1
        elif hf >= heteroplasmy_min:
            het += 1
        elif hf > 0:
            low += 1
    return n, homo, het, low


SUMMARY_HEADER = ["Sample", "Reads after trimming", "mtDNA reads", "mtDNA reads (%)",
                  "mtDNA coverage (%)", "Mean depth",
                  "Best predicted haplogroup(s)", "N. of variants",
                  "N. of homoplasmic variants",
                  "N. of heteroplasmic variants",
                  "N. of low-level variants",
                  "N. of prioritized variants"]


def summarize_samples(sample_inputs, prioritized_counts, out_file,
                      homoplasmy_threshold=0.97, heteroplasmy_min=0.03,
                      min_depth=5):
    """ Write one summary row per sample.

    Args:
        sample_inputs: list of dicts with keys sample, best_results,
            vcf (single-sample VCF), coverage (samtools depth -a) and
            mt_reads (written by write_mt_read_fraction)
        prioritized_counts: dict sample -> n. of prioritized variants
        out_file: output tsv path
        homoplasmy_threshold: min HF of homoplasmic variants
        heteroplasmy_min: min HF of heteroplasmic variants (below:
            low-level variants)
        min_depth: min depth for a position to count as covered
    """
    with open(out_file, "w") as out:
        out.write("# homoplasmic: HF >= {0}; heteroplasmic: {1} <= HF < {0}; "
                  "low-level: HF < {1}\n".format(homoplasmy_threshold,
                                                 heteroplasmy_min))
        out.write("\t".join(SUMMARY_HEADER) + "\n")
        for s in sample_inputs:
            breadth, mean_depth = read_coverage(s["coverage"], min_depth)
            n, homo, het, low = count_variants(s["vcf"], homoplasmy_threshold,
                                               heteroplasmy_min)
            row = [s["sample"], *read_mt_read_fraction(s["mt_reads"]), breadth, mean_depth,
                   read_best_haplogroups(s["best_results"]),
                   n, homo, het, low,
                   prioritized_counts.get(s["sample"], 0)]
            out.write("\t".join(str(x) for x in row) + "\n")
