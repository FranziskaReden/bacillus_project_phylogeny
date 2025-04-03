configfile: "config.yaml"

rule all:
    input:
        expand("bacillus_assemblies/prokka/{genome}/{genome}.gff", genome=config["genome"]),
        expand("alignments/faa/muscle/{marker}.faa.afa", marker=config["marker_genes"])

rule annotate:
    input: "bacillus_assemblies/{genome}/{genome}.fasta"
    output: "bacillus_assemblies/prokka/{genome}/{genome}.gff"
    shell:"""
        prokka --outdir bacillus_assemblies/prokka/{wildcards.genome} --prefix {wildcards.genome} {input} --force"
    """

rule get_gtdb:
    output: 
        "gtdb/bac120_marker_genes_all.tar.gz",
        "gtdb/bac120_marker_genes_reps.tar.gz",
        "bacillus_assemblies/prokka/hmm_PGAP.tsv",
        "gtdb/bac120_taxonomy.tsv",
        "gtdb/bac120.tree"
    script: "scripts/retrieve_gtdb_data.py"

rule retrieve_marker_genes:
    input: 
        "gtdb/bac120_marker_genes_all.tar.gz",
        "gtdb/bac120_marker_genes_reps.tar.gz",
        "bacillus_assemblies/prokka/hmm_PGAP.tsv",
        "gtdb/bac120_taxonomy.tsv"
    output: "gtdb/bac120_taxonomy_expanded.tsv"
    script: "scripts/retrieve_marker_genes.py"

rule compare_marker_genes:
    input: "gtdb/bac120_taxonomy_expanded.tsv"
    output: "bacillus_assemblies/prokka/{genome}/{genome}_closest_genomes.tsv",
        "bacillus_assemblies/prokka/{genome}/{genome}_genes_identity.tsv"
    shell: """
        python scripts/compare_marker_genes.py -a {wildcards.genome}
    """

rule choose_genomes:
    input: 
        expand("bacillus_assemblies/prokka/{genome}/{genome}_closest_genomes.tsv", genome=config["genome"]),
        expand("bacillus_assemblies/prokka/{genome}/{genome}_genes_identity.tsv", genome=config["genome"])
    output: 
        "marker_genes/best_hits.tsv",
        expand("alignments/faa/{marker}.faa", marker=config["marker_genes"]),
        expand("alignments/fna/{marker}.fna", marker=config["marker_genes"])
    script: "scripts/choose_genomes.py"

rule create_aa_alignment:
    input: "alignments/faa/{marker}.faa"
    output: 
        alignments="alignments/faa/muscle/{marker}.faa.efa",
        confidence="alignments/faa/muscle/{marker}.faa.cefa",
        best_alignment="alignments/faa/muscle/{marker}.faa.afa"
    threads: 2
    shell: """
        muscle5 -align {input} -stratified -output {output.alignments} -threads 2 > {output.alignments}.log 2>&1
        muscle5 -addconfseq {output.alignments} -output {output.confidence} >> {output.alignments}.log 2>&1
        muscle5 -maxcc {output.alignments} -output {output.best_alignment} >> {output.alignments}.log 2>&1
    """

rule overlay_alignment:
    input: "alignments/faa/muscle/{marker}.faa.afa"
    output: "alignments/fna/muscle/{marker}.fna.afa"
    shell: """
        python scripts/overlap_alignment.py -m {wildcards.marker}
    """

rule find_sub_model:
    input: "alignments/fna/muscle/{marker}.fna.afa"
    output: "alignments/fna/muscle/{marker}.iqtree"
    shell: """
        raxml-ng --msa alignments/fna/muscle/{wildcards.marker}.fna.afa --check -m GTR+G
        if [ -f alignments/fna/muscle/{wildcards.marker}.fna.afa.raxml.reduced.phy ]; then 
            iqtree2mod -s alignments/fna/muscle/{wildcards.marker}.fna.afa.raxml.reduced.phy -m MF --seqtype DNA --prefix {wildcards.marker} > alignments/fna/muscle/{wildcards.marker}.iqtree.log 2>&1
        else;
            iqtree2mod -s alignments/fna/muscle/{wildcards.marker}.fna.afa -m MF --seqtype DNA --prefix {wildcards.marker} > alignments/fna/muscle/{wildcards.marker}.iqtree.log 2>&1
        fi;
    """

rule ml_tree:
    input: 
        model_file="alignments/fna/muscle/{marker}.iqtree"
    output: "alignments/fna/muscle/{marker}.mlsearch"
    threads: 8 
    shell: """
        model=$(grep "Best-fit model according to BIC:" {input.model_file} | awk '{print $NF}')
        if [ -f alignments/fna/muscle/{wildcards.marker}.fna.afa.raxml.reduced.phy ]; then 
            ali_file=$"alignments/fna/muscle/{wildcards.marker}.fna.afa.raxml.reduced.phy"
        else
            ali_file=$"alignments/fna/muscle/{wildcards.marker}.fna.afa"
        fi;
        raxml-ng --msa $ali_file --model $model --prefix alignments/fna/muscle/{wildcards.marker}.mlsearch --all  > alignments/fna/muscle/{wildcards.marker}.mlsearch.log 2>&1
    """