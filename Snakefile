# Snakemake workflow for Bacillus project
# Copyright (c) 2025 Franziska Reden
# This workflow is licensed under the GNU General Public License v3.0 (GPL-3.0)
# See LICENSE file for details.

# Config file holds the names of the assemblies as well as marker genes.
configfile: "config.yaml"

rule all:
    """
    Will create the species tree.
    """
    input:
        "species_tree/species_tree.nwk"

rule annotate:
    """
    In the first step, the assembies as annotated using prokka.
    """
    input: "bacillus_assemblies/{genome}/{genome}.fasta"
    output: "bacillus_assemblies/prokka/{genome}/{genome}.gff"
    conda: "environments/prokka_envs.yaml"
    shell:"""
        prokka --outdir bacillus_assemblies/prokka/{wildcards.genome} --prefix {wildcards.genome} {input} --force"
    """

rule get_gtdb:
    """
    This rule retrieves all relevant GTDB data i.e. marker genes of the representative and all genomes, 
    taxonomy file, tree.
    """
    output: 
        "gtdb/bac120_marker_genes_all.tar.gz",
        "gtdb/bac120_marker_genes_reps.tar.gz",
        "bacillus_assemblies/prokka/hmm_PGAP.tsv",
        "gtdb/bac120_taxonomy.tsv",
        "gtdb/bac120.tree"
    conda: "environments/prokka_envs.yaml"
    script: "scripts/retrieve_gtdb_data.py"

rule retrieve_marker_genes:
    """
    In this step, the relevant marker genes (of the genera of interest) are retreived from the GTDB data.
    Additionally, the marker genes of the annotated assemblies are fetched.
    """
    input: 
        "gtdb/bac120_marker_genes_all.tar.gz",
        "gtdb/bac120_marker_genes_reps.tar.gz",
        "bacillus_assemblies/prokka/hmm_PGAP.tsv", 
        "gtdb/bac120_taxonomy.tsv",
        "config_genera.csv"
    output: expand("marker_genes/reps/faa/{marker}.faa", marker=config["marker_genes"]),
        expand("marker_genes/reps/fna/{marker}.fna", marker=config["marker_genes"]),
        expand("marker_genes/all/faa/{marker}.faa", marker=config["marker_genes"]),
        expand("marker_genes/all/fna/{marker}.fna", marker=config["marker_genes"])
    conda: "environments/prokka_envs.yaml"
    script: "scripts/retrieve_marker_genes.py"

rule identity_alignment:
    """
    Align aa sequences using muscle5 to calculate the average indenity between genomes based on the marker genes.
    """
    input:
        faa_ali="marker_genes/reps/faa/{marker_genes}.faa",
    output:
        faa_alignment="marker_genes/reps/faa/alignment/{marker_genes}.aln"
    threads: 1
    shell:"""
        sed -i 's/*//g' {input.faa_ali}
        bin/muscle5 -align {input.faa_ali} -output {output.faa_alignment} -threads 1
    """

rule overlay_identity_alignment:
    """
    Overlay the protein alignment onto the dna sequences to create a codon aware alignment.
    """
    input:
        faa_ali = "marker_genes/reps/faa/alignment/{marker_genes}.aln",
        fna_ali = "marker_genes/reps/fna/{marker_genes}.fna"
    output:
        fna_alignment = "marker_genes/reps/fna/alignment/{marker_genes}.aln"
    conda: "environments/prokka_envs.yaml"
    shell:"""
        python scripts/overlay_alignment.py --faa {input.faa_ali} --fna {input.fna_ali} --out {output.fna_alignment}
    """

rule choose_genomes:
    """
    In this step, the genomes to be used in the downstream analysis are chosen. Genomes are clustered 
    according to their pairwise identity (calculated on the created alignment). Genomes with an idenity of 
    over 95% are clustered and all but one representative retained for further analysis.
    Furthermore, genomes that are close to the assemblies (also based on sequence identity), but not neccessarly
    representative genomes are also included.
    """
    input:
        expand("marker_genes/reps/fna/alignment/{marker_genes}.aln", marker_genes=config["marker_genes"])
    output:
        expand("alignments/faa/{marker}.faa", marker=config["marker_genes"]),
        expand("alignments/fna/{marker}.fna", marker=config["marker_genes"])
    conda: "environments/prokka_envs.yaml"
    shell: """
        python scripts/choose_genomes_new.py
    """

rule alignment:
    """
    In this rule the chosen genomes are aligned using muscle5.
    """
    input:
        "alignments/faa/{marker}.faa"
    output:
        alignment="alignments/faa/{marker}.faa.efa",
        confidence="alignments/faa/{marker}.faa.cefa",
        best_alignment="alignments/faa/{marker}.faa.afa"
    conda: "environments/prokka_envs.yaml"
    threads: 2
    shell: """
        sed -i 's/*//g' {input}

        muscle5 -align {input} -stratified -output {output.alignment} -threads 2 > {output.alignment}.log 2>&1
        muscle5 -addconfseq {output.alignments} -output {output.confidence} >> {output.alignments}.log 2>&1
        muscle5 -maxcc {output.alignments} -output {output.best_alignment} >> {output.alignments}.log 2>&1
    """

rule overlay_alignment:
    """
    Overlay the protein alignment onto the dna sequences to create a codon aware alignment.
    """
    input:
        faa_ali = "alignments/faa/{marker}.faa.afa",
        fna_ali = "alignments/fna/{marker}.fna"
    output:
        fna_alignment = "alignments/fna/{marker_genes}.fna.aln"
    conda: "environments/prokka_envs.yaml"
    shell:"""
        python scripts/overlay_alignment.py --faa {input.faa_ali} --fna {input.fna_ali} --out {output.fna_alignment}
    """ 

rule find_sub_model:
    """
    In this step, the best fitting model is selected for each marker gene.
    """
    input: "alignments/fna/{marker}.fna.aln"
    output: "model/{marker}/{marker}.model"
    threads: 8
    shell: """
        mkdir -p model
        mkdir -p model/{wildcards.marker}

        bin/iqtree -s {input} \
                    -m MF -T AUTO \
                    --seqtype DNA \
                    --prefix model/{wildcards.marker}/{wildcards.marker}

        MODEL=$(grep "Best-fit model according to BIC:" model/{wildcards.marker}/{wildcards.marker}.iqtree | awk '{print $NF}')
        less $MODEL > {output}
    """

rule infer_ml_trees:
    """
    Infer one ML tree for each marker gene and replicate (1–50).
    """
    input: 
        model="model/{marker}/{marker}.model",
        alignment="alignments/fna/{marker}.fna.aln"
    output: "gene_trees/{marker}/{marker}_{rep}.treefile"
    threads: 10
    shell: """
        MODEL=$(< {input.model})

        bin/iqtree2 -s {input.alignment} \
                    -m "$MODEL" \
                    -T AUTO \
                    --seqtype DNA \
                    --keep-ident \
                    --prefix gene_trees/{wildcards.marker}/{wildcards.marker}_{wildcards.rep} \
        >> gene_trees/{wildcards.marker}/{wildcards.marker}_{wildcards.rep}.iqlog 2>&1
    """

rule au_test:
    """
    Run AU and other statistical test to see which trees in the 50 ML tree set are relevant.
    """
    input:
        alignment="alignments/fna/{marker}.fna.aln",
        model="model/{marker}/{marker}.model",
        trees=expand("gene_trees/{{marker}}/{{marker}}_{i}.treefile", i=range(1,51))
    output:
        trees="au_test/{marker}/{marker}.nwk",
        iqtree="au_test/{marker}/{marker}_au_test.iqtree"
    threads: 10
    shell:"""
        MODEL=$(< {input.model})

        if [ -f au_test/{wildcards.marker}/{wildcards.marker}.nwk ]; then 
            rm au_test/{wildcards.marker}/{wildcards.marker}.nwk	
        fi;

        for i in gene_trees/{wildcards.marker}/{wildcards.marker}*.treefile; do
            less $i >> au_test/{wildcards.marker}/{wildcards.marker}.nwk	
        done;

        bin/iqtree2 -s {input.alignment} \
                    -z au_test/{wildcards.marker}.nwk \
                    -m "$MODEL" \
                    -n 0 -zb 10000 -zw -au \
                    --keep-ident \
                    --prefix au_test/{wildcards.marker}/{wildcards.marker}_au_test \
                    >> au_test/{wildcards.marker}/{wildcards.marker}_au_test.aulog 2>&1
        """

rule extract_trees:
    """
    Extract the ML trees that passed at least one signifigance test.
    """
    input: 
        iqtree="au_test/{marker}/{marker}_au_test.iqtree",
        trees="au_test/{marker}/{marker}.nwk"
    output: "au_test/{marker}/{marker}_selected.nwk"
    conda: "environments/prokka_envs.yaml"
    shell:
        """
        python scripts/extract_trees.py --iqtree {input.iqtree} --trees {input.trees} --ouput {output}
        """
    
rule infer_species_tree:
    """
    Finally, create a species tree based on all significant ML trees using ASTRAL.
    """
    input:
        expand("au_test/{marker}/{marker}_selected.nwk", marker=config["marker_genes"])
    output: "species_tree/species_tree.nwk"
    threads: 16
    shell:"""
        mkdir -p species_tree

        if [ -f species_tree/all_trees.nwk ]; then
            rm species_tree/all_trees.nwk
        fi;

        for i in au_test/*/*_selected.nwk; do
            less $i > species_tree/all_trees.nwk
        done;

        bin/ASTER-Linux/bin/astral4 -i species_tree/all_trees.nwk \
                                -o species_tree/species_tree \
                                --support 2 -r 16 -s 16 -t 16 \
                                --genelength 1300 \
                                > species_tree/species_tree.log 2>&1
        """

rule rename_species_tree:
    """
    Rename the taxa in the species tree and extract support values
    """
    input: 
        treefile="species_tree/species_tree.nwk",
        taxonomy="gtdb/bac120_taxonomy_expanded.tsv"
    output: "species_tree/species_tree_support_rename.nwk"
    conda: "environments/prokka_envs.yaml"
    shell:"""
        python scripts/rename_tree.py
        """

rule get_ballistic_distances:
    input: 
        tree="species_tree/species_tree.nwk"
    output: "species_tree/ballistic_distances.tsv"
    conda: "environments/prokka_envs.yaml"
    script: "scripts/get_ballistic_distances.py"