#!/usr/bin/env python
""" Compare mtDNA variants between samples of the same individual.

Optional post-processing of MToolBox results, not part of the pipeline.

Subcommands
-----------
numt-alleles
    List the mtDNA alleles carried by NUMTs: extract the NUMT sequences
    (BED of NUMT coordinates + unmasked nuclear genome), align them to the
    mtDNA reference with GMAP and record every difference.

compare
    For every candidate variant (any allele in the pipeline's merged VCF),
    count reference and alternative reads in every sample directly from
    its alignments, so that evidence below the calling threshold is seen
    too, and report per sample: depth, strand counts, heteroplasmy (HF)
    with a 95% Wilson confidence interval, a call that accounts for depth
    (the lower confidence bound must reach --min-hf and the allele must be
    above the background noise estimated from the other individuals), the
    lowest HF the sample could detect at that position, and whether the
    allele is carried by NUMTs or recurs in other individuals. Writes a
    long table (one row per sample and variant), a per-individual table
    (one row per variant) and per-individual counts of sharing patterns
    (e.g. tumour-only, tumour+cfDNA). Also checks that the samples of an
    individual agree: same (or nested) best haplogroup, and no variant
    homoplasmic in one sample missing from another (a sample swap or mix-up).

Samples are grouped by individual using --samples (a tsv with columns
sample, individual and type) or, by default, by splitting sample names at
the last underscore (P01_tumour -> individual P01, type tumour).
"""
import argparse
import csv
import math
import os
import re
import subprocess
import sys
from collections import OrderedDict, defaultdict

import numpy as np
from scipy.stats import binom

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
from modules.prioritization import read_vcf  # noqa: E402

Z = 1.959964  # 95% confidence
ANNOTATIONS = ("Locus", "AaChange", "NtVarH", "Pathogenicity", "Clinvar")
STATUS_ORDER = ("homoplasmic", "present", "present_strand_bias", "trace",
                "absent", "not_informative", "not_called")


# ---------------------------------------------------------------- statistics

def wilson_interval(k, n, z=Z):
    """ 95% Wilson score interval of a proportion k/n. """
    if n == 0:
        return 0.0, 1.0
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def min_alt_reads(depth, min_hf):
    """ Smallest number of alternative reads whose Wilson lower bound
    reaches min_hf at this depth, or None if none does. """
    for k in range(depth + 1):
        if wilson_interval(k, depth)[0] >= min_hf:
            return k
    return None


def min_detectable_hf(depth, min_hf, power=0.8):
    """ Lowest true HF that would be called (lower confidence bound >=
    min_hf) with the given probability at this depth; 1.0 if none. """
    if depth == 0:
        return 1.0
    k = min_alt_reads(depth, min_hf)
    if k is None:
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo + hi) / 2
        if binom.sf(k - 1, depth, mid) >= power:
            hi = mid
        else:
            lo = mid
    return hi


def benjamini_hochberg(pvalues):
    """ Benjamini-Hochberg adjusted p-values, in the input order. """
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.minimum(ranked, 1.0)
    return adjusted


# ------------------------------------------------------------------- samples

def read_sample_sheet(path, samples):
    """ sample -> (individual, type), from a tsv or from sample names. """
    groups = OrderedDict()
    if path:
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                groups[row["sample"]] = (row["individual"], row["type"])
        missing = [s for s in samples if s not in groups]
        if missing:
            sys.exit("Samples missing from {}: {}".format(path, ", ".join(missing)))
        return OrderedDict((s, groups[s]) for s in samples)
    for s in samples:
        individual, _, sample_type = s.rpartition("_")
        groups[s] = (individual or s, sample_type)
    return groups


# -------------------------------------------------------------- read counts

BASES = "ACGT"


def genome_wide_error(counts, min_hf):
    """ Sequencing error rate per alternative base: fraction of bases
    differing from the majority one, over positions without any allele at
    HF >= min_hf (i.e. without real variants), divided by 3. """
    total = counts["+"] + counts["-"]
    depth = total.sum(axis=0)
    major = total.max(axis=0)
    minor_max = np.sort(total, axis=0)[-2]
    ok = (depth > 0) & (minor_max < min_hf * np.maximum(depth, 1))
    if not ok.any():
        return 0.0
    return float((depth[ok] - major[ok]).sum() / depth[ok].sum() / 3)


def count_bases(bam_file, contig, min_base_quality):
    """ Per-position base counts on each strand: returns {'+': array(4, L),
    '-': array(4, L)} with rows in ACGT order (0-based positions). """
    import pysam
    counts = {}
    with pysam.AlignmentFile(bam_file) as bam:
        if contig not in bam.references:
            sys.exit("{} not found in {}".format(contig, bam_file))
        for strand, reverse in (("+", False), ("-", True)):
            def keep(read, reverse=reverse):
                return (not (read.is_unmapped or read.is_secondary or read.is_qcfail
                             or read.is_duplicate) and read.is_reverse == reverse)
            counts[strand] = np.array(bam.count_coverage(
                contig, quality_threshold=min_base_quality, read_callback=keep))
    return counts


# ------------------------------------------------------------ NUMT alleles

def read_fasta(path):
    seqs, name = OrderedDict(), None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.upper())
    return OrderedDict((k, "".join(v)) for k, v in seqs.items())


def extract_numts(bed, genome, out_fasta, samtools="samtools"):
    """ Extract the NUMT sequences in bed (0-based) from genome. """
    if not os.path.exists(genome + ".fai"):
        subprocess.run([samtools, "faidx", genome], check=True)
    with open(genome + ".fai") as fai:
        names = {l.split("\t")[0] for l in fai}
    regions = []
    with open(bed) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) < 3 or line.startswith(("#", "track", "browser")):
                continue
            chrom = fields[0]
            if chrom not in names and chrom.startswith("chr") and chrom[3:] in names:
                chrom = chrom[3:]          # BED in UCSC names, genome in Ensembl names
            if chrom not in names:
                continue
            regions.append("{}:{}-{}".format(chrom, int(fields[1]) + 1, fields[2]))
    if not regions:
        sys.exit("None of the BED regions match sequence names in " + genome)
    with open(out_fasta, "w") as out:
        for i in range(0, len(regions), 500):
            out.write(subprocess.run([samtools, "faidx", genome] + regions[i:i + 500],
                                     check=True, stdout=subprocess.PIPE,
                                     universal_newlines=True).stdout)
    return len(regions)


def sam_differences(sam_lines, mt_seq):
    """ Differences of aligned sequences from the mtDNA reference.

    Returns {(pos, ref, alt): set of query names}, positions 1-based;
    substitutions only (ref and alt one base each).
    """
    L = len(mt_seq)
    diffs = defaultdict(set)
    for line in sam_lines:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10 or int(f[1]) & 4 or f[9] == "*":
            continue
        ref_pos, read_pos, seq = int(f[3]) - 1, 0, f[9].upper()
        for length, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
            length = int(length)
            if op in "M=X":
                for i in range(length):
                    r = mt_seq[(ref_pos + i) % L]
                    q = seq[read_pos + i]
                    if q != r and q in BASES and r in BASES:
                        diffs[((ref_pos + i) % L + 1, r, q)].add(f[0])
                ref_pos += length
                read_pos += length
            elif op in "IS":
                read_pos += length
            elif op in "DN":
                ref_pos += length
    return diffs


def numt_alleles(args):
    mt = list(read_fasta(args.mt_fasta).values())[0]
    numt_fasta = args.numt_fasta
    if numt_fasta is None:
        if not (args.bed and args.genome):
            sys.exit("Give --numt-fasta, or --bed and --genome (unmasked).")
        numt_fasta = os.path.splitext(args.out)[0] + "_sequences.fasta"
        n = extract_numts(args.bed, args.genome, numt_fasta)
        print("Extracted {} NUMT sequences to {}".format(n, numt_fasta))
    gmap = subprocess.run(["gmap", "-D", args.gmap_db_dir, "-d", args.gmap_db, "-f", "samse",
                           "-n", "1", "-t", str(args.threads), numt_fasta],
                          check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                          universal_newlines=True)
    diffs = sam_differences(gmap.stdout.splitlines(), mt)
    with open(args.out, "w") as out:
        out.write("pos\tref\talt\tn_numts\tnumts\n")
        for (pos, ref, alt), names in sorted(diffs.items()):
            out.write("{}\t{}\t{}\t{}\t{}\n".format(pos, ref, alt, len(names),
                                                    ",".join(sorted(names))))
    print("{} NUMT alleles written to {}".format(len(diffs), args.out))


def load_numt_alleles(path):
    if not path:
        return {}
    with open(path) as fh:
        return {(int(r["pos"]), r["alt"]): int(r["n_numts"])
                for r in csv.DictReader(fh, delimiter="\t")}


# -------------------------------------------------------------- concordance

def read_haplogroups(path):
    """ sample -> list of best haplogroups, from the pipeline's
    results/haplogroups/<mt>_<n>_best_results.csv; {} if not found. """
    if not path or not os.path.exists(path):
        return {}
    haplogroups = {}
    with open(path) as fh:
        for row in csv.reader(fh):
            if row and row[0] != "SampleID":
                haplogroups[row[0]] = [h for h in (row[1] if len(row) > 1 else "").split(";") if h]
    return haplogroups


def haplogroup_relation(a, b):
    """ 'identical', 'nested' (one is a sub-haplogroup of the other by
    name, e.g. H1 and H1a), 'different' or 'unknown', for two lists of
    best haplogroups. """
    if not a or not b:
        return "unknown"
    if set(a) == set(b):
        return "identical"
    if any(x.startswith(y) or y.startswith(x) for x in a for y in b):
        return "nested"
    return "different"


def homoplasmic_discordance(hf_depth, min_depth=20, homoplasmy=0.97, absent_below=0.5):
    """ Homoplasmic variants on which a sample disagrees with the other
    samples of the same individual.

    For each variant homoplasmic (HF >= homoplasmy) in at least one sample,
    the samples with at least min_depth reads are split into homoplasmic and
    absent (HF < absent_below); the minority group disagrees with the
    individual (both groups when they are the same size, e.g. two samples).

    Args:
        hf_depth: {variant: {sample: (hf, depth)}} for one individual
    Returns:
        {sample: (n_compared, [variants on which it disagrees])}
    """
    samples = sorted({s for by_sample in hf_depth.values() for s in by_sample})
    result = {s: [0, []] for s in samples}
    for variant, by_sample in hf_depth.items():
        covered = {s: hf for s, (hf, depth) in by_sample.items() if depth >= min_depth}
        homoplasmic = [s for s, hf in covered.items() if hf >= homoplasmy]
        if not homoplasmic:
            continue
        absent = [s for s, hf in covered.items() if hf < absent_below]
        for s in covered:
            result[s][0] += 1
        if not absent:
            continue
        if len(homoplasmic) > len(absent):
            minority = absent
        elif len(absent) > len(homoplasmic):
            minority = homoplasmic
        else:
            minority = homoplasmic + absent
        for s in minority:
            result[s][1].append(variant)
    return {s: (n, d) for s, (n, d) in result.items()}


def write_concordance(prefix, samples, groups, haplogroups, snv_rows, args):
    """ Per-sample haplogroup and homoplasmic-variant concordance with the
    other samples of the same individual. Returns the samples to check. """
    by_individual = defaultdict(list)
    for s in samples:
        by_individual[groups[s][0]].append(s)
    to_check = []
    with open(prefix + "_concordance.tsv", "w") as out:
        out.write("individual\tsample\ttype\thaplogroup\tindividual_haplogroup\t"
                  "haplogroup_relation\thomoplasmic_compared\thomoplasmic_discordant\t"
                  "discordant_positions\tverdict\n")
        for ind, members in by_individual.items():
            # the individual's haplogroup: the most frequent among its samples
            calls = [tuple(haplogroups.get(s, [])) for s in members]
            known = [c for c in calls if c]
            reference = list(max(set(known), key=known.count)) if known else []
            hf_depth = {}
            for v, ev in snv_rows:
                hf_depth["{}{}>{}".format(v["pos"], v["ref"], v["alt"])] = {
                    s: (ev[s]["hf"], ev[s]["depth"]) for s in members}
            discordance = homoplasmic_discordance(hf_depth, args.concordance_min_depth,
                                                  args.homoplasmy)
            for s in members:
                hap = haplogroups.get(s, [])
                relation = haplogroup_relation(hap, reference) if len(members) > 1 else "single_sample"
                n, discordant = discordance.get(s, (0, []))
                if len(discordant) > args.max_discordant:
                    verdict = "likely_mismatch"
                elif relation == "different":
                    verdict = "check_haplogroup"
                else:
                    verdict = "ok"
                if verdict != "ok":
                    to_check.append((s, verdict))
                out.write("\t".join([ind, s, groups[s][1], ";".join(hap), ";".join(reference),
                                     relation, str(n), str(len(discordant)),
                                     ",".join(discordant[:20]), verdict]) + "\n")
    return to_check


# ------------------------------------------------------------------ compare

def vcf_contig(vcf):
    """ Sequence name of the first record of a (possibly gzipped) VCF. """
    import gzip
    with open(vcf, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    with (gzip.open(vcf, "rt") if gz else open(vcf)) as fh:
        for line in fh:
            if not line.startswith("#"):
                return line.split("\t")[0]
    return None


def candidate_variants(vcf):
    """ Alleles of the merged VCF: list of dicts with pos, ref, alt, kind,
    annotations and per-sample VCF calls. """
    _, samples, records = read_vcf(vcf)
    variants = []
    for rec in records:
        kind = "snv" if len(rec["ref"]) == 1 and len(rec["alt"]) == 1 else "indel"
        if rec["alt"] in ("*", "."):
            continue
        variants.append({"pos": rec["pos"], "ref": rec["ref"], "alt": rec["alt"],
                         "kind": kind, "calls": rec["calls"],
                         "annotations": {k: rec["info"].get(k, "") for k in ANNOTATIONS}})
    return samples, variants


def _float(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def compare(args):
    samples, variants = candidate_variants(args.vcf)
    unknown = set(args.exclude) - set(samples)
    if unknown:
        sys.exit("--exclude: not in the VCF: " + ", ".join(sorted(unknown)))
    samples = [s for s in samples if s not in set(args.exclude)]
    groups = read_sample_sheet(args.samples, samples)
    numts = load_numt_alleles(args.numt_alleles)
    vcf_base = re.sub(r"(\.annotated)?\.vcf(\.gz)?$", "", os.path.basename(args.vcf))
    pattern = args.bam_pattern or ("results/{sample}/map/{sample}_" + vcf_base +
                                   "_OUT-sorted.realign.bam")
    snvs = [v for v in variants if v["kind"] == "snv"]
    contig = args.contig or vcf_contig(args.vcf)

    # read counts at every candidate SNV, in every sample
    evidence = {}  # (sample, i) -> dict
    sample_error = {}  # sample -> sequencing error rate per alternative base
    for s in samples:
        bam = pattern.format(sample=s)
        if not os.path.exists(bam):
            sys.exit("Alignment not found: {} (use --bam-pattern)".format(bam))
        counts = count_bases(bam, contig, args.min_base_quality)
        sample_error[s] = genome_wide_error(counts, args.min_hf)
        for i, v in enumerate(snvs):
            p = v["pos"] - 1
            a = BASES.index(v["alt"]) if v["alt"] in BASES else None
            depth = int(counts["+"][:, p].sum() + counts["-"][:, p].sum())
            fwd = int(counts["+"][a, p]) if a is not None else 0
            rev = int(counts["-"][a, p]) if a is not None else 0
            evidence[(s, i)] = {"depth": depth, "alt_fwd": fwd, "alt_rev": rev}
        print("Counted reads in {} (error rate per base: {:.2g})".format(s, sample_error[s]),
              file=sys.stderr)

    # first pass: confidence intervals and clear presence, for the background
    for (s, i), e in evidence.items():
        alt = e["alt_fwd"] + e["alt_rev"]
        e["alt"] = alt
        e["hf"] = alt / e["depth"] if e["depth"] else 0.0
        e["ci_low"], e["ci_high"] = wilson_interval(alt, e["depth"])

    # background noise of each allele: pooled over the samples of the other
    # individuals where it is not clearly present
    individuals = sorted({g[0] for g in groups.values()})
    for i in range(len(snvs)):
        for ind in individuals:
            alt = depth = 0
            for s in samples:
                e = evidence[(s, i)]
                if groups[s][0] != ind and e["ci_low"] < args.min_hf:
                    alt += e["alt"]
                    depth += e["depth"]
            err = (alt / depth) if depth else 0.0
            for s in samples:
                if groups[s][0] == ind:
                    # the sample's own error rate is a floor: the other
                    # individuals can show no error read at a position by chance
                    evidence[(s, i)]["background"] = max(err, sample_error[s], args.min_error)

    # significance above background, BH-adjusted within each sample
    for s in samples:
        pvals = []
        for i in range(len(snvs)):
            e = evidence[(s, i)]
            e["p"] = float(binom.sf(e["alt"] - 1, e["depth"], e["background"])) if e["alt"] else 1.0
            pvals.append(e["p"])
        for i, q in enumerate(benjamini_hochberg(pvals)):
            evidence[(s, i)]["q"] = float(q)

    mdhf_cache = {}
    def mdhf(depth):
        if depth not in mdhf_cache:
            mdhf_cache[depth] = min_detectable_hf(depth, args.min_hf)
        return mdhf_cache[depth]

    for (s, i), e in evidence.items():
        e["min_detectable_hf"] = mdhf(e["depth"])
        significant = e["q"] < args.alpha
        if e["hf"] >= args.homoplasmy and e["ci_low"] >= args.min_hf:
            e["status"] = "homoplasmic"
        elif e["ci_low"] >= args.min_hf and significant:
            both = e["alt_fwd"] >= args.min_strand_reads and e["alt_rev"] >= args.min_strand_reads
            e["status"] = "present" if both else "present_strand_bias"
        elif significant:
            e["status"] = "trace"
        else:
            e["status"] = "absent"

    # indels: values of the pipeline VCF
    indels = [v for v in variants if v["kind"] == "indel"]
    indel_evidence = {}
    for j, v in enumerate(indels):
        for s in samples:
            call = v["calls"].get(s, {})
            hf, dp, low = _float(call.get("HF")), _float(call.get("DP")), _float(call.get("CILOW"))
            e = {"depth": int(dp) if dp else 0, "alt_fwd": "", "alt_rev": "", "alt": "",
                 "hf": hf if hf is not None else 0.0, "ci_low": low if low is not None else 0.0,
                 "ci_high": _float(call.get("CIUP"), ""), "background": "", "q": "",
                 "min_detectable_hf": mdhf(int(dp)) if dp else 1.0}
            if hf is None:
                e["status"] = "not_called"
            elif hf >= args.homoplasmy:
                e["status"] = "homoplasmic"
            elif e["ci_low"] >= args.min_hf:
                e["status"] = "present"
            else:
                e["status"] = "trace"
            indel_evidence[(s, j)] = e

    all_rows = [(v, {s: evidence[(s, i)] for s in samples}) for i, v in enumerate(snvs)]
    all_rows += [(v, {s: indel_evidence[(s, j)] for s in samples}) for j, v in enumerate(indels)]

    # recurrence: individuals in which each allele is present
    present_in = {}
    for v, ev in all_rows:
        present_in[(v["pos"], v["alt"])] = {groups[s][0] for s, e in ev.items()
                                            if e["status"] in ("present", "homoplasmic")}

    write_outputs(args.out, samples, groups, all_rows, present_in, numts, args)

    haplogroup_file = args.haplogroups or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(args.vcf))), "haplogroups",
        vcf_base + "_best_results.csv")
    to_check = write_concordance(args.out, samples, groups, read_haplogroups(haplogroup_file),
                                 [(v, {s: evidence[(s, i)] for s in samples})
                                  for i, v in enumerate(snvs)], args)
    print("Wrote {}_concordance.tsv".format(args.out))
    for s, verdict in to_check:
        print("WARNING: {} {}: check {}_concordance.tsv".format(s, verdict, args.out),
              file=sys.stderr)


def write_outputs(prefix, samples, groups, rows, present_in, numts, args):
    out_dir = os.path.dirname(prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    called = {"present", "present_strand_bias", "homoplasmic"}
    types_of = defaultdict(list)
    for s in samples:
        types_of[groups[s][0]].append((s, groups[s][1]))

    long_cols = ["individual", "sample", "type", "pos", "ref", "alt", "kind", "depth",
                 "alt_fwd", "alt_rev", "hf", "ci_low", "ci_high", "background_error",
                 "q_value", "min_detectable_hf", "status", "numt_alleles",
                 "other_individuals"] + list(ANNOTATIONS)
    wide_rows, overlap = [], defaultdict(lambda: defaultdict(int))
    with open(prefix + "_long.tsv", "w") as out:
        out.write("\t".join(long_cols) + "\n")
        for v, ev in rows:
            key = (v["pos"], v["alt"])
            for ind, members in types_of.items():
                calls = {s: ev[s] for s, _ in members}
                if not any(e["status"] in called for e in calls.values()):
                    continue
                max_hf = max(e["hf"] for e in calls.values() if e["status"] in called)
                others = len(present_in[key] - {ind})
                for s, sample_type in members:
                    e = calls[s]
                    status = e["status"]
                    if status == "absent" and e["min_detectable_hf"] > max_hf:
                        status = "not_informative"
                    e["final_status"] = status
                    out.write("\t".join(_fmt(x) for x in [
                        ind, s, sample_type, v["pos"], v["ref"], v["alt"], v["kind"],
                        e["depth"], e["alt_fwd"], e["alt_rev"], e["hf"], e["ci_low"],
                        e["ci_high"], e["background"], e["q"], e["min_detectable_hf"],
                        status, numts.get(key, 0), others]
                        + [v["annotations"][k] for k in ANNOTATIONS]) + "\n")
                statuses = [calls[s]["final_status"] for s, _ in members]
                if len(members) == 1:
                    # nothing to compare with: the pattern of e.g. a cfDNA-only
                    # individual must not count as a cfDNA_only variant
                    sharing = "single_sample"
                elif all(st == "homoplasmic" for st in statuses):
                    sharing = "homoplasmic_all"
                else:
                    sharing = "+".join(sorted({t for (s, t) in members
                                               if calls[s]["final_status"] in called})) + "_only"
                    if len(members) > 1 and all(calls[s]["final_status"] in called for s, _ in members):
                        sharing = "all"
                flags = []
                if numts.get(key):
                    flags.append("numt_allele")
                if others >= args.recurrent:
                    flags.append("recurrent")
                if any(calls[s]["final_status"] == "present_strand_bias" for s, _ in members):
                    flags.append("strand_bias")
                if any(calls[s]["final_status"] == "not_informative" for s, _ in members):
                    flags.append("low_power")
                overlap[ind][sharing] += 1
                wide_rows.append((ind, v, members, calls, sharing, flags, others))

    with open(prefix + "_by_individual.tsv", "w") as out:
        types = sorted({t for g in groups.values() for t in [g[1]]})
        out.write("\t".join(["individual", "pos", "ref", "alt", "kind", "sharing", "flags",
                             "other_individuals"] +
                            ["{}_{}".format(t, c) for t in types for c in ("hf", "depth", "status")] +
                            list(ANNOTATIONS)) + "\n")
        for ind, v, members, calls, sharing, flags, others in wide_rows:
            by_type = {t: calls[s] for s, t in members}
            cells = []
            for t in types:
                e = by_type.get(t)
                cells += [_fmt(e["hf"]), _fmt(e["depth"]), e["final_status"]] if e else ["", "", ""]
            out.write("\t".join([ind, str(v["pos"]), v["ref"], v["alt"], v["kind"], sharing,
                                 ",".join(flags), str(others)] + cells +
                                [v["annotations"][k] for k in ANNOTATIONS]) + "\n")

    with open(prefix + "_overlap.tsv", "w") as out:
        out.write("individual\tsharing\tn_variants\n")
        for ind in sorted(overlap):
            for sharing, n in sorted(overlap[ind].items()):
                out.write("{}\t{}\t{}\n".format(ind, sharing, n))
    print("Wrote {0}_long.tsv, {0}_by_individual.tsv and {0}_overlap.tsv".format(prefix))


def _fmt(x):
    if isinstance(x, float):
        return "{:.4g}".format(x)
    return str(x)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command")

    n = sub.add_parser("numt-alleles", help="list the mtDNA alleles carried by NUMTs")
    n.add_argument("--bed", help="NUMT coordinates (BED, 0-based)")
    n.add_argument("--genome", help="UNMASKED nuclear genome fasta (indexed or indexable)")
    n.add_argument("--numt-fasta", help="NUMT sequences, instead of --bed and --genome")
    n.add_argument("--mt-fasta", required=True, help="mtDNA reference used by the pipeline")
    n.add_argument("--gmap-db-dir", default="gmap_db", help="default: gmap_db")
    n.add_argument("--gmap-db", required=True, help="GMAP database of the mtDNA reference, "
                   "e.g. rCRS (built by the pipeline in gmap_db/)")
    n.add_argument("--threads", type=int, default=4)
    n.add_argument("--out", default="numt_alleles.tsv")

    c = sub.add_parser("compare", help="compare variants between samples of an individual")
    c.add_argument("--vcf", required=True,
                   help="merged (annotated) VCF, e.g. results/vcf/rCRS_GRCh38.annotated.vcf")
    c.add_argument("--samples", help="tsv with columns sample, individual, type "
                   "(default: split sample names at the last underscore)")
    c.add_argument("--bam-pattern", help="alignments, with {sample} (default: "
                   "results/{sample}/map/{sample}_<mt>_<n>_OUT-sorted.realign.bam)")
    c.add_argument("--contig", help="mtDNA sequence name in the alignments (default: the VCF's)")
    c.add_argument("--numt-alleles", help="output of numt-alleles")
    c.add_argument("--min-hf", type=float, default=0.03,
                   help="a call needs the lower 95%% confidence bound of HF >= this (default 0.03)")
    c.add_argument("--homoplasmy", type=float, default=0.97, help="default 0.97")
    c.add_argument("--alpha", type=float, default=0.01,
                   help="FDR for alleles above background noise (default 0.01)")
    c.add_argument("--min-strand-reads", type=int, default=2,
                   help="alternative reads needed on each strand (default 2)")
    c.add_argument("--min-base-quality", type=int, default=25, help="default 25")
    c.add_argument("--min-error", type=float, default=1e-4,
                   help="lowest background error rate allowed (default 1e-4)")
    c.add_argument("--recurrent", type=int, default=2,
                   help="flag alleles present in at least this many other individuals")
    c.add_argument("--exclude", nargs="+", default=[], metavar="SAMPLE",
                   help="samples to leave out, e.g. ones flagged in _concordance.tsv")
    c.add_argument("--haplogroups", help="best haplogroups per sample (default: "
                   "results/haplogroups/<mt>_<n>_best_results.csv next to the VCF's folder)")
    c.add_argument("--concordance-min-depth", type=int, default=20,
                   help="depth needed to compare a homoplasmic variant between samples (default 20)")
    c.add_argument("--max-discordant", type=int, default=2,
                   help="homoplasmic variants a sample may lack before it is flagged as a "
                   "likely mismatch with its individual (default 2)")
    c.add_argument("--out", default="comparison/mt_comparison", help="output prefix")

    args = parser.parse_args(argv)
    if args.command == "numt-alleles":
        numt_alleles(args)
    elif args.command == "compare":
        compare(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
