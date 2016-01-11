deBWT:
BWT string parallel construction method based on de Bruijn graph theory.

Introduction:
deBWT is an innovative BWT string construction method based on graph theory, whose most parts can be implemented in parallel. Therefore, this method is extremely fast, especially for high similar sequence. For example, it can build ten human genomes in 2 hours under 32 threads (Exclude two txt transfer parts that need about 2 hours. This part can be eliminated as long as we can handle jellyfish's binary file). Based on our big data tests, deBWT is faster at least several folds than any other methods around the world.
deBWT is open source and free for non-commercial use.
deBWT is mainly designed by Bo Liu & Dixian Zhu, optimized and developed by Dixian Zhu in the Center for Bioinformatics, Harbin Institute of Technology, China.

Memory usage:
In order to keep k-mer counting in high performance, we use 120G RAM to do this part. And based on our data, 30G human genome sequence, which has about 5G distinct k-mers, the software at most need 120G RAM. Please make sure your disk has enough space for temporary files (will be deleted after running). The requirement based on your data. For our 30G data, we need around 200G disk to store those files.

Getting start:
Download the source code from https://github.com/DixianZhu/deBWT or https://github.com/hitbc/deBWT
cd deBWT; make
./deBWT -o target_bwt_file -j jellyfish_directory -t thread_number sequence_file 

Usage:
./deBWT [options] reference
options:
-t: thread number (default: 32)
-o: bwt string output object file
-j: your jellyfish directory
reference: sequence in fasta or fastq format
