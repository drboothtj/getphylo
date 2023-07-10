# getphylo: GEnbank To PHYLOgeny
a python package for automated generation of heuristic phylogenetic trees from genbank files

## Description
`getphylo` was designed to automatically build multi-locus phylogenetic trees from GenBank files. The workflow consists of the following steps: i) extract protein coding sequences; ii) screen for suilable markers; iii) align individual marker sequences and create a combined alignment; and iv) produce a tree from the combined alignment. Please see the 'parameters' section below for a full list of parameters.

## Installation

**Important:** The following intructuons are for Linux and assume you have Python and pip installed already. For installation on different systems, please see the links provided for this individual dependences.

### A. Manual Installation

#### 1. Install getphylo
##### Installation with PyPI
The easiest way to install `getphylo` is using the command: 

`pip install getphylo` 

This will fetch and install the latest version form: https://pypi.org/project/getphylo/

#### Install getphylo's dependencies
#### 2. Installation of dependences

**Important:** Getphylo requires DIAMOND, MUSCLE and FastTree2 to be installed to work correctly. These must be installed manually. BioPython should be installed automatically. Below are instructions to install each dependency:

##### 2. Install DIAMOND

You can install DIAMOND with the following command:

`sudo apt install diamond-aligner`

Further instructions for installing DIAMOND can be found, here: https://github.com/bbuchfink/diamond/wiki/2.-Installation.


##### 3. Install MUSCLE

You can install MUSCLE with the following command:

`sudo apt install muscle`

Further instructions for installing MUSCLE can be found, here: http://www.drive5.com/muscle/muscle.html

**IMPORTANT:** getphylo is not currently compatable with MUSCLE 5. We recommend using version 3.8 as this was the version used in development. Although, any version below 5.0 should work without issues.


##### 4. Install FastTree2

You can install FastTree2 with the following command:

`sudo apt install fasttree`

Further instructions for installing FastTree2 can be found, here: http://www.microbesonline.org/fasttree/.


### B. Conda Installation











#### Installation with Conda
#####


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
