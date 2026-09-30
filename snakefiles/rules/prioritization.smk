# Variant prioritization and per-sample summary (human only).
# Needs human_analysis_tab from rules/haplogroup.smk and the annotated VCF
# from rules/annotation.smk.
from modules.prioritization import prioritize_variants, summarize_samples

def human_samples(ref_genome_mt, ref_genome_n):
    tab = human_analysis_tab
    return list(tab.loc[(tab["ref_genome_mt"] == ref_genome_mt) &
                        (tab["ref_genome_n"] == ref_genome_n), "sample"])

def prioritization_inputs(wildcards):
    samples = human_samples(wildcards.ref_genome_mt, wildcards.ref_genome_n)
    refs = {"ref_genome_mt": wildcards.ref_genome_mt, "ref_genome_n": wildcards.ref_genome_n}
    base = "results/{sample}/haplogroup/{sample}_{ref_genome_mt}_{ref_genome_n}"
    return {
        "annotated_vcf": "results/vcf/{ref_genome_mt}_{ref_genome_n}.annotated.vcf".format(**refs),
        "merged_diffs": [(base + "_merged_diff.csv").format(sample=s, **refs) for s in samples],
        "best_results": [(base + "_best_results.csv").format(sample=s, **refs) for s in samples],
        "vcfs": ["results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.vcf.gz".format(sample=s, **refs)
                 for s in samples],
        "coverage": ["results/{sample}/map/{sample}_{ref_genome_mt}_{ref_genome_n}_OUT-sorted.realign.bam.cov".format(
                     sample=s, **refs) for s in samples]
    }

rule prioritize:
    input:
        unpack(prioritization_inputs)
    output:
        prioritized = "results/prioritization/{ref_genome_mt}_{ref_genome_n}_prioritized_variants.txt",
        summary = "results/prioritization/{ref_genome_mt}_{ref_genome_n}_summary.txt"
    params:
        samples = lambda wildcards: human_samples(wildcards.ref_genome_mt, wildcards.ref_genome_n),
        homoplasmy_threshold = config.get("prioritization", {}).get("homoplasmy_threshold", 0.97),
        heteroplasmy_min = config.get("prioritization", {}).get("heteroplasmy_min", 0.03),
        min_depth = config["mtvcf_main_analysis"]["minrd"]
    message: "Prioritizing variants for {wildcards.ref_genome_mt}_{wildcards.ref_genome_n}"
    run:
        counts = prioritize_variants(list(zip(params.samples, input.merged_diffs)),
                                     input.annotated_vcf, output.prioritized)
        sample_inputs = [{"sample": s, "best_results": b, "vcf": v, "coverage": c}
                         for s, b, v, c in zip(params.samples, input.best_results,
                                               input.vcfs, input.coverage)]
        summarize_samples(sample_inputs, counts, output.summary,
                          homoplasmy_threshold=params.homoplasmy_threshold,
                          heteroplasmy_min=params.heteroplasmy_min,
                          min_depth=params.min_depth)
