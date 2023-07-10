# getphylo: GEnbank To PHYLOgeny
a python package for automated generation of heuristic phylogenetic trees from genbank files

## Description
getphylo was designed to automatically build multi-locus phylogenetic trees from GenBank files. The workflow consists of the following steps: i) extract protein coding sequences; ii) screen for suilable markers; iii) align individual marker sequences and create a combined alignment; and iv) produce a tree from the combined alignment. Please see the 'parameters' section below for a full list of parameters.

## Installation

**Important:** The following intructuons are for Linux and assume you have Python and pip installed already. For installation on different systems, please see the links provided for this individual dependences.

### A. Manual Installation

#### 1. Install getphylo
##### Installation with PyPI
The easiest way to install `getphylo` is using the command: 

`pip install getphylo`

This will fetch and install the latest version from: https://pypi.org/project/getphylo/

#### 2. Installation of dependences

**Important:** getphylo requires DIAMOND, MUSCLE and FastTree2 to be installed to work correctly. These must be installed manually. BioPython should be installed automatically. Below are instructions to install each dependency:

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

##### 1. Create a new conda environment from the getphylo.yml file
WIP

##### 2. Install getphylo
See part A.1 for instructions.


## Usage
This package has been designed to be as easy to run as possible. Simply navigate to a working directory containing .gbk files and input:

`getphylo`

This will run the software on default settings and output all alignments and trees.

There may be occasions where you need to change the default settings. A full list of options and flags can be viewed with:

`getphylo -h`

Below is a breif describtion of each parameter.

### Parameters
The following is a complete list of parameters that can be used when running getphylo:


#### Help and Logging
`-h, --help            show this help message and exit`
`-l {CRITICAL,FATAL,ERROR,WARN,WARNING,INFO,DEBUG,NOTSET}, --logging {CRITICAL,FATAL,ERROR,WARN,WARNING,INFO,DEBUG,NOTSET}
                        set the logging level (default: ERROR)`

EXPLAIN

#### Input and Output 
  ´-g GBKS, --gbks GBKS  string indicating the genbank files to use in the phylogeny (default: *.gbk)´
  ´-o OUTPUT, --output OUTPUT
                        a string designating the name of the folder to output the results(default: output)´
  ´-s SEED, --seed SEED  path to a genbankfile with for the target organism (default: None)´
  
  ´-r RANDOM_SEED_NUMBER, --random-seed-number RANDOM_SEED_NUMBER
                        interger to be used as a seed for randomising loci selection, random if left as None(default: None)´
  ´-t TAG, --tag TAG     string indicating the GenBank annotations to extract (default: locus_tag)´

#### Loci Thresholding Parameters
  ´-f FIND, --find FIND  integer indicating the number of loci to find in the seed genome (default: -1)´
  ´-max MAXLENGTH, --maxlength MAXLENGTH
                        interger indicating the minimum length of loci to be included in the analysis (default: 2000)´
  ´-min MINLENGTH, --minlength MINLENGTH
                        interger indicating the minimum length of loci to be included in the analysis (default: 200)´
  ´-minl MINLOCI, --minloci MINLOCI
                        minimum number of loci required to continue to alignment and tree building steps (default: 1)´
  ´-maxl MAXLOCI, --maxloci MAXLOCI
                        maximum number of loci required to continue to alignment and tree building steps (default: 1000)´
  ´-p PRESENCE, --presence PRESENCE
                        interger indicating the percentage of genomes each loci must be present in (default: 100)´

#### Tree Building
`-b BUILD_ALL, --build-all BUILD_ALL
                        build phylogenetic trees for all loci, not just concatenated alignment (default: 0)`

#### Handeling Poorly Formatted Data
`-ia IGNORE_BAD_ANNOTATIONS, --ignore-bad-annotations IGNORE_BAD_ANNOTATIONS
                        ignore missing annotations - NOT RECCOMMENDED (default: False)`
`-ir IGNORE_BAD_RECORDS, --ignore-bad-records IGNORE_BAD_RECORDS
                        ignore poorly formated records - NOT RECCOMMENDED (default: False)`

#### Checkpointing
  `-cp {START,FASTA_EXTRACTED,DIAMOND_BUILT,SINGLETONS_IDENTIFIED,SINGLETONS_SEARCHED,SINGLETONS_THRESHOLDED,SINGLETONS_EXTRACTED,SINGLETONS_ALIGNED,ALIGNMENTS_COMBINED,TREES_BUILT,DONE}, --checkpoint {START,FASTA_EXTRACTED,DIAMOND_BUILT,SINGLETONS_IDENTIFIED,SINGLETONS_SEARCHED,SINGLETONS_THRESHOLDED,SINGLETONS_EXTRACTED,SINGLETONS_ALIGNED,ALIGNMENTS_COMBINED,TREES_BUILT,DONE}
                        string indicating the checkpoint to start from START = default
                         FASTA_EXTRACTED = Skip extracting fasta sequences from genbank files
                         DIAMOND_BUILT = Skip building diamond databases
                         SINGLETONS_IDENTIFIED = Skip identifying singletons from the seed genome
                         SINGLETONS_SEARCHED = Skip searching singletons against other genomes
                         SINGLETONS_THRESHOLDED = Skip thresholding of singletons
                         SINGLETONS_EXTRACTED = Skip extract fasta sequences for alignments
                         SINGLETONS_ALIGNED = Skip individual protein alignments
                         ALIGNMENTS_COMBINED = Skip combining alignments
                         TREES_BUILT = Skip building trees
                         DONE = Done
                         (default: START)`

#### Performance                       
`-c CPUS, --cpus CPUS  The number of cpus to use for paralleslisation (default: 1)`
  
  



## Examples
See the ´example_data´ folder in this directory for three example trees produced by `getphylo`.

1. A phylogeny of bacterial genomes
2. A phylogeny of eukaryotic mitochondrial genomes
3. A phylogeny of a biosynthetic super cluster
4. A phylogeny of primate genomes
