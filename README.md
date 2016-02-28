# de Bruijn branch-based BWT constructor (deBWT)  #


----------

## Introduction: ##
- DeBWT is a novel parallel BWT construction approach. It innovatively represents and organizes the suffixes of input sequence with a novel data structure, de Bruijn branch encode. This data structure takes the advantage of de Bruijn graph to facilitate the comparison between the suffixes with long common prefix, which breaks the bottleneck of the BWT construction of repetitive genomic sequences. Meanwhile, deBWT also utilizes the structure of de Bruijn graph for reducing unnecessary comparisons between suffixes.

- DeBWT has good scalability to construct BWT in parallel computing. It is well-suited to run on multiple core servers or clusters to construct the BWT of large collections of genome sequences. 

- DeBWT needs to install and run in Linux. It is open source and free for non-commercial use.
- DeBWT is designed by Bo Liu & Dixian Zhu, developed by Dixian Zhu in the Center for Bioinformatics, Harbin Institute of Technology, China.

==============================================================================

## Memory usage: ##

- We employ Jellyfish tool (http://www.genome.umd.edu/jellyfish.html) to implement the k-mer counting step of deBWT. To obtain high performance, this step requires large RAM space. For example, it allocates 120G RAM to count the k-mers of 10 human genomes. Moreover, please also make sure your disk has enough space for temporary files (will be deleted after running). The requirement depends on your data. For a collection of 10 human genomes (about 30 Gbp in total), it needs around 200G hard disk space to store temporary files.

==============================================================================

## Getting start: 

- As deBWT calls Jellyfish to implement the k-mer counting step, Jellyfish needs to be installed before use deBWT. A tuned configuration of Jellyfish has been built in the current version of deBWT. The user does not need to tune the parameters of Jellyfish by himself.
- The source code of deBWT is available at https://github.com/DixianZhu/deBWT or https://github.com/hitbc/deBWT. 

- Compile the source code with the following commands:

		cd deBWT; make
		./deBWT -o target_bwt_file -j jellyfish_directory -t thread_number sequence_file 

==============================================================================

## Usage: ##

./deBWT [options] input_sequence
                **(make sure the sequence don't contain any uncertain characters like 'N')**
### options: ###
- -t: number of thread (default: 8)
- -o: bwt string output file
- -j: the directory of jellyfish

### input_sequence: ###

- the sequence to be indexed, should be in fasta format
