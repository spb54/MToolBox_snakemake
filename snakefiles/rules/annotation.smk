# Functional annotation of merged VCFs with mtoolnote.
# Needs reference_tab and config from variant_calling.snakefile.
from modules.config_parsers import get_analysis_species, mtoolnote_species

rule annotate_vcf:
    input:
        vcf = "results/vcf/{ref_genome_mt}_{ref_genome_n}.vcf"
    output:
        vcf = "results/vcf/{ref_genome_mt}_{ref_genome_n}.annotated.vcf"
    params:
        species = lambda wildcards: mtoolnote_species(
            get_analysis_species(wildcards.ref_genome_mt, reference_tab=reference_tab,
                                 config_species=config["species"])
        )
    log: log_dir + "/vcf/{ref_genome_mt}_{ref_genome_n}_annotate_vcf.log"
    message: "Annotating {input.vcf} with mtoolnote (species: {params.species})"
    # run in its own process: mtoolnote's sqlite connection must not be shared
    # between snakemake threads
    shell:
        """
        python -c 'import sys, mtoolnote; mtoolnote.annotate(sys.argv[1], sys.argv[2], species=sys.argv[3])' \
            {input.vcf} {output.vcf} {params.species} &> {log}
        """
