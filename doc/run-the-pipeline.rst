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
- :code:`results/prioritization/<ref_genome_mt>_<ref_genome_n>_summary.txt`: per sample coverage, mean depth, best haplogroup(s), number of variants (homoplasmic, heteroplasmic above/below :code:`prioritization: hf_threshold` in :code:`config.yaml`, default 0.8) and of prioritized variants
