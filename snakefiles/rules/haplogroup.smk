# Haplogroup prediction (human only) with mt_classifier and Phylotree r17.
# Needs analysis_tab, reference_tab, config and source_dir from
# variant_calling.snakefile.
import shutil

from modules.config_parsers import get_haplo_prediction_files, get_human_analyses

human_analysis_tab = get_human_analyses(analysis_tab, reference_tab=reference_tab,
                                        config_species=config["species"])
muscle_exe = config.get("haplogroup", {}).get("muscle_exe") or shutil.which("muscle")
classifier_data_dir = os.path.join(str(source_dir), "data", "classifier")

rule main_mt_hpred:
    input:
        single_fasta = "results/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}.fasta"
    output:
        haplo_pred = "results/{sample}/haplogroup/{sample}_{ref_genome_mt}_{ref_genome_n}.csv",
        haplo_pred_sorted = "results/{sample}/haplogroup/{sample}_{ref_genome_mt}_{ref_genome_n}.sorted.csv",
        merged_diff = "results/{sample}/haplogroup/{sample}_{ref_genome_mt}_{ref_genome_n}_merged_diff.csv",
        best_results = "results/{sample}/haplogroup/{sample}_{ref_genome_mt}_{ref_genome_n}_best_results.csv"
    params:
        basename = "results/{sample}/haplogroup/{sample}_{ref_genome_mt}_{ref_genome_n}",
        muscle_exe = muscle_exe,
        data_dir = classifier_data_dir,
        source_dir = source_dir
    log: log_dir + "/{sample}/{sample}_{ref_genome_mt}_{ref_genome_n}_main_mt_hpred.log"
    message: "Predicting haplogroup of {wildcards.sample} from {input.single_fasta}"
    shell:
        """
        PYTHONPATH={params.source_dir}${{PYTHONPATH:+:$PYTHONPATH}} python -m modules.mt_classifier \
            -i {input.single_fasta} \
            -m {params.muscle_exe} \
            -b {params.basename} \
            -s {output.best_results} \
            -d {params.data_dir} \
            -n {wildcards.sample} &> {log}
        """

rule merge_haplo_best_results:
    input:
        lambda wildcards: get_haplo_prediction_files(human_analysis_tab,
                                                     ref_genome_mt=wildcards.ref_genome_mt,
                                                     ref_genome_n=wildcards.ref_genome_n,
                                                     suffix="_best_results.csv")
    output:
        best_results = "results/haplogroups/{ref_genome_mt}_{ref_genome_n}_best_results.csv"
    message: "Collecting best haplogroup predictions for {wildcards.ref_genome_mt}"
    run:
        with open(output.best_results, "w") as out:
            out.write("SampleID,Best predicted haplogroup(s)\n")
            for best_results in input:
                with open(best_results) as fh:
                    out.write(fh.read())
