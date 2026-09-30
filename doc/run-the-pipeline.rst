Run MToolBox
============

MToolBox is made by several snakemake workflows which can be run independently. We provide wrappers for the most common tasks and analyses. Using these wrappers will save a lot of typing and headache for the lazy users (probably *you*). Cool, isn't it? :)
 
All the wrappers accepts snakemake arguments and parse automatically the `config.yaml` configuration file required by snakemake.

Setting up a working directory
------------------------------

Replace :code:`/path/to/MToolBox/dir/` with the MToolBox installation path and :code:`/path/to/analysis/dir` with the folder where you wish to run your analysis.

.. code-block:: bash
    
    export MTOOLBOX_DIR=/path/to/MToolBox/dir/
    
    cd /path/to/analysis/dir
    cp $MTOOLBOX_DIR/config.yaml .
    cp $MTOOLBOX_DIR/cluster.yaml .
    # create folders needed by the workflow
    mkdir -p data/reads
    mkdir -p data/genomes
    mkdir -p logs/cluster_jobs

Copy read datasets to `data/reads` and genome fasta files to `data/genomes`.

.. code-block:: bash
    
    # some code

How to run the MToolBox wrappers
--------------------------------

Running the wrappers is as simple as this:

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-<wrapper> <snakemake arguments>

*E.g.* if you want to run the MToolBox-variant-calling wrapper and print the commands it will execute, you can run

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-variant-calling-bash -p

You can also display a graphical representation of the workflow by running

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-variant-calling-bash --dag | display

This will show the workflow in a browser. Alternatively, you can save the workflow representation in a file by running

.. code-block:: bash
    
    export PATH=/path/to/MToolBox/dir/:$PATH
    
    MToolBox-variant-calling-bash --dag > workflow.svg

Available wrappers
------------------

- `MToolBox-variant-calling`_
- `MToolBox-variant-annotation`_
- `MToolBox-human-haplogroup-prediction`_
- `MToolBox-human-prioritization`_

Each wrapper runs variant calling first, so any of them can be started from raw reads.

MToolBox-variant-calling
^^^^^^^^^^^^^^^^^^^^^^^^

Performs QC, quality trimming of raw reads, read alignment, alignment filtering, variant calling. The final output is a VCF file.

In the VCF, :code:`HF` is the heteroplasmy fraction of each alternative allele, i.e. its read count divided by :code:`DP`. Alternative alleles are only counted from bases with quality >= :code:`mtvcf_main_analysis: Q` and at least :code:`tail_mismatch` bases from the read ends (to avoid end-of-read artefacts), and :code:`DP` counts the reads passing the same filters at that position, so that :code:`HF` is not biased by the read-end filter. For indels, only informative reads are counted in both :code:`DP` and the allele count: reads covering the base preceding the indel and the first base after the repeat or homopolymer containing it (e.g. the poly-C tract at 303-315 or the CA repeat at 514-523), since reads ending inside the repeat cannot show whether the indel is present. Indel alleles that no informative read supports are dropped.

MToolBox-variant-annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotates the merged VCF of each reference genome pair with `mtoolnote <https://github.com/mitoNGS/mtoolnote>`_, writing :code:`results/vcf/<ref_genome_mt>_<ref_genome_n>.annotated.vcf`.
The species passed to mtoolnote is the :code:`species` column of :code:`data/reference_genomes.tab`, or :code:`species` in :code:`config.yaml` if set. Use :code:`human` (:code:`hsapiens` is also accepted) for human, or one of the non-human species supported by mtoolnote (e.g. :code:`ggallus`, :code:`mmusculus`); non-human annotation queries Ensembl BioMart, so it needs internet access.

MToolBox-human-haplogroup-prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Human only (analyses whose reference species is human). Assigns a Phylotree (build 17) haplogroup to each sample's consensus fasta with mt-classifier, which aligns it to RSRS with MUSCLE 3.8. Outputs, per sample, in :code:`results/<sample>/haplogroup/`:

- :code:`<sample>_<ref_genome_mt>_<ref_genome_n>.csv` and :code:`.sorted.csv`: haplogroup predictions
- :code:`..._merged_diff.csv`: variants relative to RSRS, the MHCS of the predicted haplogroup and rCRS
- :code:`..._best_results.csv`: best predicted haplogroup(s)

and :code:`results/haplogroups/<ref_genome_mt>_<ref_genome_n>_best_results.csv` with the best haplogroup(s) of all samples.

MToolBox-human-prioritization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Human only. Runs the complete analysis: variant calling, haplogroup prediction, annotation, then

- :code:`results/prioritization/<ref_genome_mt>_<ref_genome_n>_prioritized_variants.txt`: variants private to a sample, i.e. differing from RSRS, from the MHCS of its haplogroup and from rCRS, with their mtoolnote annotations, sorted by ascending nucleotide variability (:code:`NtVarH`)
- :code:`results/<sample>/<sample>_<ref_genome_mt>_<ref_genome_n>_mtDNA_reads.tsv`: number of reads after trimming (all libraries of the sample, both mates of surviving pairs plus surviving unpaired reads), number of reads mapped to the mtDNA (primary alignments, before duplicate removal, so that duplicates are counted in both) and their percentage
- :code:`results/prioritization/<ref_genome_mt>_<ref_genome_n>_summary.txt`: per sample reads after trimming, mtDNA reads and their percentage, coverage, mean depth, best haplogroup(s), number of variants, split by heteroplasmy fraction (HF) into homoplasmic (HF >= :code:`prioritization: homoplasmy_threshold` in :code:`config.yaml`, default 0.97), heteroplasmic (HF >= :code:`heteroplasmy_min`, default 0.03) and low-level, and number of prioritized variants
