# Bacillota Phylogenomics Workflow

This repository contains a Snakemake workflow to reconstruct a species tree for *Bacillota* genomes using marker genes of the GTDB and maximum likelihood inference using IQ-Tree2. The workflow integrates annotation, marker gene extraction, alignment, model selection, gene tree inference, AU testing, and species tree reconstruction.

## Workflow Overview

1. **Annotate genomes** with **Prokka** ([Seemann, 2014](https://doi.org/10.1093/bioinformatics/btu153)).  
2. **Retrieve GTDB data** ([Parks et al., 2020](https://doi.org/10.1038/s41587-020-0501-8)) and extract marker genes.  
3. **Align sequences** with **MUSCLE v5** ([Edgar, 2022](https://doi.org/10.1093/molbev/msac025)) and filter genomes >95% identical.  
4. **Infer gene trees** with **IQ-TREE2** ([Minh et al., 2020](https://doi.org/10.1093/molbev/msaa015)) and select supported trees via AU tests.  
5. **Build species tree** using **ASTRAL** ([Mirarab & Warnow, 2015](https://doi.org/10.1093/bioinformatics/btv234)).  
6. **Post-process tree**: rename taxa, extract support values, and compute ballistic distances.  

## Installation

This workflow requires [Snakemake](https://snakemake.readthedocs.io/) and [Conda / Mamba](https://mamba.readthedocs.io/) (recommended).

Clone the repository:

```bash
git clone https://github.com/FranziskaReden/bacillus_project_phylogeny.git
cd bacillus_project_phylogeny
```

Additionally, the following third-party software must be installed manually (or downloaded as binaries) and placed in the `bin/` directory of this repository: [**IQ-TREE2**](http://www.iqtree.org/), [**MUSCLE5**](https://github.com/rcedgar/muscle), [**ASTRAL**](https://github.com/smirarab/ASTRAL). After download, ensure the binaries are executable and available in your bin:

```bash
chmod +x bin/iqtree2
chmod +x bin/muscle5
chmod +x bin/ASTER-Linux/bin/astral4
```

The workflow was tested using the following versions: **Prokka v1.14.6** ([Seemann, 2014](https://doi.org/10.1093/bioinformatics/btu153)), **MUSCLE v5.1** ([Edgar, 2022](https://doi.org/10.1093/molbev/msac025)), **IQ-TREE2 v2.4.0** ([Minh et al., 2020](https://doi.org/10.1093/molbev/msaa015)), **ASTRAL IV** ([Mirarab & Warnow, 2015](https://doi.org/10.1093/bioinformatics/btv234)).

## Run the workflow

The assemblies to be studied have to be placed into the [bacillus_assemblies](bacillus_assemblies) folder. The names of the marker genes and assemblies have to be stated in the [config.yaml](config.yaml) file. Also, the [config_genera.csv](config_genera.csv) holds the closest putative relatives of the assemblies.

Run snakemake with the following command: 

```bash
snakemake --use-conda --cores <N>
```

Note, that the `--use-conda` command is necessary as the workflow uses a conda environment as specified in [environments](environments).

## Output

- **Species tree:** `species_tree/species_tree.nwk`
- **AU test trees:** `au_test/{marker}/{marker}_au_test.iqtree`
- **Alignments:** `alignments/faa/` and `alignments/fna/`
- **Additional intermediate files** are stored in the respective folders: `marker_genes/`, `model/`, `gene_trees/`

## License

This workflow is licensed under the **Creative Commons Licence (CC0-1.0 license)**. See the [LICENSE](LICENSE) file for details.

## References

1. Edgar, R. C. (2022). MUSCLE v5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. *Nature Communications, 13*, 6968. [https://doi.org/10.1038/s41467-022-34630-w](https://doi.org/10.1038/s41467-022-34630-w)

2. Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. *Molecular Biology and Evolution, 37*(5), 1530–1534. [https://doi.org/10.1093/molbev/msaa015](https://doi.org/10.1093/molbev/msaa015)

3. Mirarab, S., & Warnow, T. (2015). ASTRAL: Genome-scale coalescent-based species tree estimation. *Bioinformatics, 30*(17), i541–i548. [https://doi.org/10.1093/bioinformatics/btv234](https://doi.org/10.1093/bioinformatics/btv234)

4. Seemann, T. (2014). Prokka: Rapid prokaryotic genome annotation. *Bioinformatics, 30*(14), 2068–2069. [https://doi.org/10.1093/bioinformatics/btu153](https://doi.org/10.1093/bioinformatics/btu153)

5. Parks, D. H., Chuvochina, M., Chaumeil, P. A., Rinke, C., Mussig, A. J., Hugenholtz, P., & Rappé, M. S. (2020). GTDB: An ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank-normalized and complete genome-based taxonomy. *Nucleic Acids Research, 50*(D1), D785–D794. [https://doi.org/10.1093/nar/gkz1019](https://doi.org/10.1093/nar/gkz1019)
