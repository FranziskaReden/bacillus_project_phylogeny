import urllib.request
import os

def retrieve_marker_genes():

    print(f"Retrieving the HMM database for PGAP...")
    # Retrieve the HMM database for PGAP
    file_name = "bacillus_assemblies/prokka/hmm_PGAP.tsv"
    if os.path.exists("bacillus_assemblies/prokka") == False:
        os.mkdir("bacillus_assemblies/prokka")
    if os.path.exists(file_name) == False:
        # Retrieve the HMM database for
        url = "ftp://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.tsv"
        urllib.request.urlretrieve(url, file_name)

    print(f"Retrieving the GTDB taxonomy...")
    file_name = "gtdb/bac120_taxonomy.tsv.gz"
    # Retrieve the taxobomy of the GTDB genomes
    if os.path.exists("gtdb") == False:
        os.mkdir("gtdb")
    if os.path.exists(file_name.replace('.gz', '')) == False:
        url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120_taxonomy.tsv.gz"
        urllib.request.urlretrieve(url, file_name)
        os.system("gunzip gtdb/bac120_taxonomy.tsv.gz")

    print(f"Retrieving the GTDB representative marker genes...")
    file_name = "gtdb/bac120_marker_genes_reps.tar.gz"
    # Retrieve marker genes of representative genomes
    if os.path.exists(file_name) == False:
        url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_reps/bac120_marker_genes_reps.tar.gz"
        urllib.request.urlretrieve(url, file_name)
        # Unzip files
        os.system(f"tar -xvf {file_name} -C gtdb")

    print(f"Retrieving all GTDB marker genes...")
    file_name = "gtdb/bac120_marker_genes_all.tar.gz"
    # Retrieve marker genes of all genomes
    if os.path.exists(file_name) == False:
        url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/genomic_files_all/bac120_marker_genes_all.tar.gz"
        urllib.request.urlretrieve(url, file_name)
        # Unzip files
        os.system(f"tar -xvf {file_name} -C gtdb")

    print(f"Retrieving the GTDB tree...")
    file_name = "gtdb/bac120.tree.gz"
    # Retrieve marker genes of representative genomes
    if os.path.exists(file_name.replace('.gz', '')) == False:
        url = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/bac120.tree.gz"
        urllib.request.urlretrieve(url, file_name)
        # Unzip files
        os.system(f"gunzip {file_name}")
    
if __name__ == "__main__":
    retrieve_marker_genes()