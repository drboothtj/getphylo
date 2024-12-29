# getphylo: GEnbank To PHYLOgeny
a python package for automated generation of heuristic phylogenetic trees from genbank files

## Description
getphylo was designed to automatically build multi-locus phylogenetic trees from GenBank files. The workflow consists of the following steps: i) extract protein coding sequences; ii) screen for suitable markers; iii) align individual marker sequences and create a combined alignment; and iv) produce a tree from the combined alignment. Please see the 'parameters' section below for a full list of parameters.

## Installation

The easiest way to install `getphylo` is using the command: 

`pip install getphylo`

This will fetch and install the latest version from: https://pypi.org/project/getphylo/

For full installation instructions, please see the [getphylo wiki](https://github.com/drboothtj/getphylo/wiki/Installation).

**Important:** getphylo requires DIAMOND(>=2.0.14.152), MUSCLE(>=3.8.1551) and FastTree2(>=2.1.11) to be installed to work correctly. These must be installed manually. Further instructions are [availiable on the wiki](https://github.com/drboothtj/getphylo/wiki/Installation).

## Quick-start
This package has been designed to be as easy to run as possible. Simply navigate to a working directory containing .gbk files and input:

`getphylo`

This will run the software with default settings.

A full list of options and flags can be viewed with:

`getphylo -h`

A [full list of parameters](https://github.com/drboothtj/getphylo/wiki/Parameter-List) and [further usage examples](https://github.com/drboothtj/getphylo/wiki/Advanced-Usage-(Case-Studies)) are availiable on the wiki.

## Example Analysis and Datasets
Example outputs and benchmarking data can be found in the [getphylo benchmarking repository](https://github.com/drboothtj/getphylo_benchmarking). The example data includes:

1. A phylogeny of bacterial genomes,
2. A phylogeny of a biosynthetic gene cluster,
3. A phylogeny of primate genomes,
4. A phylogeny of Eurotiomycete fungi.

## Citation
If you use `getphylo`, please cite:

> Booth, T. J., Shaw, S., & Weber, T. (2023). getphylo: rapid and automatic generation of multi-locus phylogenetic trees. BioRxiv, 2023.07.26.550493. 

DOI: https://doi.org/10.1101/2023.07.26.550493

## Patch Notes
### Version  0
- 0.1.0 
	- beta version initial release
- 0.1.1 
	- added support for MUSCLE5
- 0.1.2 
	- now raises an error if translations are present but empty
	- error messages from the extract module are now more informative
	- fixed a fatal issue with --build-all
- 0.2.0
	- now supports iqtree using the --method parameter
- 0.2.1
	- now able to provide custom paths for binary dependencies
	- parser now has argument groups and is more readable
	- file exists error message more informative
- 0.2.2
  - added error message when users attempt to input directory instead of a search string
- 0.3.0
  - now supports modifying blastp thresholds, including parameters for identity and coverage
  - fixed typos in parser
  - fixed crashing when provided with directories with spaces in the names
- 0.3.1
  - fixed issue with the query and subject cover in diamond
- 0.3.2
  - added version info to setup.py and README for none-python dependencies
  -
