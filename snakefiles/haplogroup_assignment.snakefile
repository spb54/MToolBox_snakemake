# Haplogroup prediction for human samples: variant calling, then
# mt_classifier on each sample's consensus fasta.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.snakefile).resolve().parent.parent))

from modules.config_parsers import get_human_analyses, get_ref_pair_files, parse_config_tabs

configfile: "config.yaml"

# the target rule must be the first rule, so it is defined before the includes
_analysis_tab, _reference_tab, _ = parse_config_tabs(analysis_tab_file="data/analysis.tab",
                                                     reference_tab_file="data/reference_genomes.tab",
                                                     datasets_tab_file="data/datasets.tab")
_human_analysis_tab = get_human_analyses(_analysis_tab, reference_tab=_reference_tab,
                                         config_species=config["species"])
if _human_analysis_tab.empty:
    sys.exit("Haplogroup prediction is only available for human: no analysis in "
             "data/analysis.tab has a human mt reference genome.")

rule all_haplo_prediction:
    input:
        get_ref_pair_files(_human_analysis_tab,
                           "{results}/haplogroups/{ref_genome_mt}_{ref_genome_n}_best_results.csv"),
        get_ref_pair_files(_analysis_tab, "{results}/summary/{ref_genome_mt}_{ref_genome_n}_summary.txt")

include: "variant_calling.snakefile"
include: "rules/haplogroup.smk"
