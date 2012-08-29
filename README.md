=========================================================================================
Tangram 0.1.0        Release Distribution Documentation               2012-08-29
Jiantao Wu (jiantaowu.xining@gmail.com), Marth Lab [1], Boston College Biology Department
=========================================================================================


Introduction
=========================

Tangram is a C/C++ command line toolbox for structural variation(SV) detection based on 
MOSAIK [2] alignments. It takes advantage of both read-pair and split-read algorithms 
and is extremely fast and memory-efficient. Powered by the Bamtools API [3], Tangram can 
call SV events on multiple BAM files (a population) simutaneously to increase the 
sensitivity on low-coverage dataset.Currently it reports mobile element insertions (MEI). 
More other SV event types will be introduced soon.


Obtaining and Compiling
=========================

> git clone git://github.com/jiantao/Tangram.git
> cd src
> make


Detection pipeline
=========================

Currently, Tangram contains four sub-programs:

1. tangram_scan   : Scan through the bam file and calculate the fragment length distribution
                    for each library in that bam file. It will output the fragment length
                    distribution files for each input bam file.

2. tangram_merge  : If more than one bam files need to be scanned, this program will combine
                    all the fragment length distribution files together. It will output the
                    merged fragment length distribution file that enable the detection of
                    multiple bam files simutaneously. This step is optional if only one bam
                    file was used.

3. tangram_index  : Index the normal and special (MEI sequences) reference file. It will output
                    the indexed refrence file. This step is required for split read algorithm.

4. tangram_detect : Detect the SV events from the MOSAIK aligned BAM files. It will output the
                    unfiltered VCF files.


The overall detection pipeline for Tangram looks like the following

tangram_scan (input: BAM files) --> fragment length distribution files  \
                                                                         \
                                                                          -----> tangram_detect (BAM files) --> VCF files
                                                                         /
tangram_index (input: reference fasta files) --> indexed reference file /

For the detailed usage of each program, please run "$PROGRAM -help"


Bug Report
=========================

Please report bugs using the built-in bug reporting feature in github or by sending the author 
an email.


References
=========================

[1] http://bioinformatics.bc.edu/marthlab/Main_Page 
[2] https://github.com/wanpinglee/MOSAIK 
[3] https://github.com/pezmaster31/bamtools
