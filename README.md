# getphylo
a python package for automated generation of heuristic phylogenetic trees from genbank files

## Description
WIP

## Installation

**Important:** The following intructuons are for Linux and assume you have Python and pip installed already. For installation on different systems, please see the links provided for this individual dependences.

### 1. Installation of getphylo

#### Installation with PyPI
Coming soon!

#### Manual Installation

If you want to install `getphylo` manually, you can do so by cloning this repository.

`git clone https://github.com/DrBoothTJ/getphylo`
`cd getphylo`
`python setup.py install`

### 2. Installation of dependences

**Important:** Getphylo requires DIAMOND, MUSCLE and FastTree2 to be installed to work correctly. These must be installed manually. BioPython should be installed automatically. Below are instructions to install each dependency:

#### Installing DIAMOND

You can install DIAMOND with the following command:

`sudo apt install diamond-aligner`

Further instructions for installing DIAMOND can be found, here: https://github.com/bbuchfink/diamond/wiki/2.-Installation.


#### Installing MUSCLE

You can install MUSCLE with the following command:

WIP


#### Installing FastTree2

You can install FastTree2 with the following command:

`sudo apt install fasttree`

Further instructions for installing FastTree2 can be found, here: http://www.microbesonline.org/fasttree/.

#### Installing BioPython

If you installed `getphylo` using PyPI, BioPython should be installed automatically. If you installed `getphylo` manually, you will need to install BioPython yourself.
This can be done using the command:

`pip install biopython`

Further instructions for installing FastTree2 can be found, here: https://biopython.org/wiki/Download.


## Usage
This package has been designed to be as easy to run as possible. Simply navigate to a working directory containing .gbk files and input:

`getphylo`

This will run the software on default settings and output all alignments and trees.

There may be occasions where you need to change the default settings. A full list of options and flags can be viewed with:

`getphylo -h`

Below is a breif describtion of each flag.

### Flags
WIP

## Examples
See the ´example_data´ folder in this directory for three example trees produced by `getphylo`.

1. A phylogeny of bacterial genomes
2. A phylogeny of eukaryotic mitochondrial genomes
3. A phylogeny of a biosynthetic super cluster
4. A phylogeny of primate genomes
