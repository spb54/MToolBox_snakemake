# Complete human analysis: variant calling, haplogroup prediction,
# functional annotation, variant prioritization and per-sample summary.
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
    sys.exit("Prioritization is only available for human: no analysis in "
             "data/analysis.tab has a human mt reference genome.")

rule all_prioritization:
    input:
        get_ref_pair_files(_human_analysis_tab,
                           "{results}/haplogroups/{ref_genome_mt}_{ref_genome_n}_best_results.csv"),
        get_ref_pair_files(_human_analysis_tab,
                           "{results}/prioritization/{ref_genome_mt}_{ref_genome_n}_prioritized_variants.txt"),
        get_ref_pair_files(_human_analysis_tab,
                           "{results}/prioritization/{ref_genome_mt}_{ref_genome_n}_summary.txt")

include: "variant_calling.snakefile"
include: "rules/haplogroup.smk"
include: "rules/annotation.smk"
include: "rules/prioritization.smk"
