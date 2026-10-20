# Per-sample summary for any species: reads, mtDNA reads and their
# percentage, coverage and variant counts. Included by variant_calling.snakefile.
# The human prioritization summary (rules/prioritization.smk) adds the
# haplogroup and prioritized variant columns.
from modules.prioritization import summarize_samples, write_mt_read_fraction

def sample_libraries(sample):
    return sorted(datasets_tab.loc[datasets_tab["sample"] == sample, "library"].unique())

rule mt_read_fraction:
    input:
        # trimmed reads: only so that this runs after trimming, the read
        # counts come from the Trimmomatic logs
        trimmed = lambda wildcards: expand("data/reads_filtered/{sample}_{library}_qc_R1.fastq.gz",
                                           sample=wildcards.sample,
                                           library=sample_libraries(wildcards.sample)),
        # mtDNA alignments before duplicate removal, as the total includes duplicates
        mt_bams = lambda wildcards: expand("results/{sample}/map/OUT_{sample}_{library}_{ref_genome_mt}_{ref_genome_n}/"
                                           "{sample}_{library}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.bam",
                                           sample=wildcards.sample,
                                           library=sample_libraries(wildcards.sample),
                                           ref_genome_mt=wildcards.ref_genome_mt,
                                           ref_genome_n=wildcards.ref_genome_n)
    output:
        mt_reads = "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}_mtDNA_reads.tsv"
    params:
        logs = lambda wildcards: expand(log_dir + "/trimmomatic/{sample}_{library}_trimmomatic.log",
                                        sample=wildcards.sample,
                                        library=sample_libraries(wildcards.sample))
    message: "Counting mtDNA reads of {wildcards.sample}"
    run:
        write_mt_read_fraction(wildcards.sample, params.logs, input.mt_bams, output.mt_reads)

def analysis_samples(ref_genome_mt, ref_genome_n):
    return list(analysis_tab.loc[(analysis_tab["ref_genome_mt"] == ref_genome_mt) &
                                 (analysis_tab["ref_genome_n"] == ref_genome_n), "sample"])

def summary_inputs(wildcards):
    samples = analysis_samples(wildcards.ref_genome_mt, wildcards.ref_genome_n)
    refs = {"ref_genome_mt": wildcards.ref_genome_mt, "ref_genome_n": wildcards.ref_genome_n}
    return {
        "vcfs": ["results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz".format(sample=s, **refs)
                 for s in samples],
        "coverage": ["results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.bam.cov".format(
                     sample=s, **refs) for s in samples],
        "mt_reads": ["results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}_mtDNA_reads.tsv".format(
                     sample=s, **refs) for s in samples]
    }

rule summary:
    input:
        unpack(summary_inputs)
    output:
        summary = "results/summary/{ref_genome_mt}_{ref_genome_n}_summary.txt"
    params:
        samples = lambda wildcards: analysis_samples(wildcards.ref_genome_mt, wildcards.ref_genome_n),
        homoplasmy_threshold = config.get("prioritization", {}).get("homoplasmy_threshold", 0.97),
        heteroplasmy_min = config.get("prioritization", {}).get("heteroplasmy_min", 0.03),
        min_depth = config["mtvcf_main_analysis"]["minrd"]
    message: "Summarizing samples of {wildcards.ref_genome_mt}_{wildcards.ref_genome_n}"
    run:
        sample_inputs = [{"sample": s, "vcf": v, "coverage": c, "mt_reads": m}
                         for s, v, c, m in zip(params.samples, input.vcfs,
                                               input.coverage, input.mt_reads)]
        summarize_samples(sample_inputs, None, output.summary,
                          homoplasmy_threshold=params.homoplasmy_threshold,
                          heteroplasmy_min=params.heteroplasmy_min,
                          min_depth=params.min_depth)
