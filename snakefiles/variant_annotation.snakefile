# Functional annotation: variant calling, then mtoolnote on the merged VCFs.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(workflow.snakefile).resolve().parent.parent))

from modules.config_parsers import get_genome_vcf_files, parse_config_tabs

configfile: "config.yaml"

# the target rule must be the first rule, so it is defined before the includes
_analysis_tab, _, _ = parse_config_tabs(analysis_tab_file="data/analysis.tab",
                                        reference_tab_file="data/reference_genomes.tab",
                                        datasets_tab_file="data/datasets.tab")

rule all_annotation:
    input:
        get_genome_vcf_files(_analysis_tab, annotation=True)

include: "variant_calling.snakefile"
include: "rules/annotation.smk"
