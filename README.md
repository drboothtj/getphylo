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

[[MORE EXAMPLES OF BASIC USAGE]]

## Parameters
The following is a complete list of parameters that can be used when running getphylo:


### Help and Logging
`-h, --help show this help message and exit`

Provides a list of parameters and other information.

`-l {CRITICAL,FATAL,ERROR,WARN,WARNING,INFO,DEBUG,NOTSET}, --logging {CRITICAL,FATAL,ERROR,WARN,WARNING,INFO,DEBUG,NOTSET}
                        set the logging level (default: ERROR)`

Set the logging level. Most users will only ever need the default setting (ERROR) or to set `-l INFO` for basic logging. `DEBUG` can also be used, but this is extremely verbose and is intended for development 
purposes only.

### Input and Output 
  `-g GBKS, --gbks GBKS  string indicating the genbank files to use in the phylogeny (default: *.gbk)`

This will point getphylo to the input files. Ensure that the input is formatted as a string, e.g.: `'-g ./input/*.gb'`. 
  
  `-o OUTPUT, --output OUTPUT
                        a string designating the name of the folder to output the results(default: output)`

This points getphylo to a specified output folder. By default results will be stored at `./output`. getphylo will make a new directory if necissary.

  `-s SEED, --seed SEED  path to a genbankfile with for the target organism (default: None)`

This sets the 'seed genome.' By default getphyl choses the first genome alphabetically in the input file list. The seed genome will only impact the analysis in two specific scenarios. Firstly, if you are 
analysing genomes of dramatically different sizes. In this case it is optimal to choose the smallest genome as this will make the run slightly faster as the starting list of possible marker genes will be lower. 
Secondly, if you have the `--presence` parameter is set lower than 100. In this case, the seed genome may effect the resulting tree. In this case, it is advisable to use the outgroup as the seed. However, in the case it
is valuable to run the analysis more than once using different seeds to ensure there is no effect. See the section on the `--presence` parameter below for more details.
  
  `-r RANDOM_SEED_NUMBER, --random-seed-number RANDOM_SEED_NUMBER
                        interger to be used as a seed for randomising loci selection, random if left as None(default: None)`
                        
If a limit is set on the number of loci wiht `--maxloci`, getphylo will randomly select marker genes. The random seed paramater can be used to choose a custom seed for this shuffling event.
                        
  `-t TAG, --tag TAG     string indicating the GenBank annotations to extract (default: locus_tag)`

This defines the genbank annotations for the protein sequences that getphylo will extract. getphylo searches for all CDS features with the provided tag. If your data does not contain `locus_tag` annotations, another 
common tag to use is `protein_id`. Ensure all data is uniformly formatted, when using getphylo!

#### Loci Thresholding Parameters
  `-f FIND, --find FIND  integer indicating the number of loci to find in the seed genome (default: -1)`
  
For large genomes, runtime can be reduced by limiting the number of loci to search for in the seed genome. This is **NOT RECCOMMENDED**.

  `-max MAXLENGTH, --maxlength MAXLENGTH
                        interger indicating the minimum length of loci to be included in the analysis (default: 2000)`

Max length can be used to limit the maximum length of marker genes. This filter is limited to exclude longer genes with multiple domains (e.g. PKSs) that may confuse the analaysis. This can be raised to 
include more loci.

  `-min MINLENGTH, --minlength MINLENGTH
                        interger indicating the minimum length of loci to be included in the analysis (default: 200)`

Max length can be used to limit the maximum length of marker genes. This filter is limited to exclude short genes (e.g. pseudogenes) with multiple domains that may confuse the analaysis. This can be lowered to 
include more loci.

  `-minl MINLOCI, --minloci MINLOCI
                        minimum number of loci required to continue to alignment and tree building steps (default: 1)`

The minimum number of marker genes required for the workflow to continue. This can be used to prematurely end runs where limited number of marker genes are avaliable. 

  `-maxl MAXLOCI, --maxloci MAXLOCI
                        maximum number of loci required to continue to alignment and tree building steps (default: 1000)`

The minimum number of marker genes required for the workflow to continue. You may wish to limit the total number to improve the performance of the analysis. For example, when analysing closely related taxa it 
is possible to find 1000s of hits; this will increase runtime significantly and may cause memory issues on some machines.
                    
  `-p PRESENCE, --presence PRESENCE
                        interger indicating the percentage of genomes each loci must be present in (default: 100)`

The percentage of genomes the marker needs to be present in. **Use with caution!** This paramater is very useful when analysing distantly related strains as there may be few markers avalible in all genomes.
Lowering the percentage has two potential drawbacks. Firstly, it will introduce missing data into the alignment which may decrease the quality of the resulting tree. Alignments should be checked by the user to assess quality. This will also make the list of markers dependent on the seed genome. This is not necisarrily a problem, but it is advisable to check the output from multiple seeds in this instance. 

#### Tree Building
`-b BUILD_ALL, --build-all BUILD_ALL
                        build phylogenetic trees for all loci, not just concatenated alignment (default: 0)`

This instructs getphylo to run fasttree on all alignments. By default a tree is built from the combined alignment only. This parameter will increase the runtime significantly.

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
  
  



## Example Analysis and Datasets
See the ´example_data´ folder in this directory for three example trees produced by `getphylo`.

1. A phylogeny of bacterial genomes
2. A phylogeny of eukaryotic mitochondrial genomes
3. A phylogeny of a biosynthetic super cluster
4. A phylogeny of primate genomes
