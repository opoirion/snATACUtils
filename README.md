# Demultiplexer for scATAC-Seq

This simple package aims to insert scATA-Seq index tags inside read ID. The barcoding and tagging strategy follows the protocol described in (Preissl et. al, 2018, doi:10.1038/s41593-018-0079-3). In Preissl et. al., Each barcode consists of four 8-bp long indexes (i5, i7, p5 and p7). The first 8 bp of Index1 correspond to the p7 barcode and the last 8 bp to the i7 barcode. The first 8 bp of Index2 correspond to the i5 barcode and the last 8 bp to the p5 barcode.

## dependencies:
    * github.com/dsnet/compress/bzip2

```bash
go get -u github.com/dsnet/compress/bzip2 # install the dependencies
```

## Installation

```bash
	git clone https://gitlab.com/Grouumf/ATACdemultiplex.git
	cd ATACdemultiplex
	go install
```

## Usage

```bash
Usage of ATACdemultiplex:
  -compressionMode int
        compressionMode for native bzip2 lib
         (1 faster -> 9 smaller) <default: 6> (default 6)
  -fastq_I1 string
        fastq index file index paired read 1
  -fastq_I2 string
        fastq index file index paired read 2
  -fastq_R1 string
        fastq read file index paired read 1
  -fastq_R2 string
        fastq read file index paired read 2
  -guess_nb_lines
        guess automatically position of the lines (for mulithread). May be not safe in some situation
  -index_no_replicate string
        <OPTIONAL> path toward indexes when only 1 replicate is used
  -index_replicate_r1 string
        <OPTIONAL> path toward indexes of R1 replicates (i.e. replicate number 1)
  -index_replicate_r2 string
        <OPTIONAL> path toward indexes of R2 replicates (i.e. replicate number 2)
  -max_nb_mistake int
        Maximum number of mistakes allowed to assign a reference read id (default 2) (default 2)
  -max_nb_reads int
        <OPTIONAL> max number of reads to process (default 0 => None)
  -nbThreads int
        number of threads to use (default 1)
  -output_tag_name string
        tag for the output file names (default None)
  -taglength int
        <OPTIONAL> number of nucleotides to consider at the end
         and begining (default 8) (default 8)
  -use_bzip2_go_lib
        use bzip2 go library instead of native C lib (slower)
  -write_extensive_logs
        write extensive logs (can consume extra RAM memory and slower the process)
  -write_logs
        write logs (might slower the execution time)

```

## multithreads

The speed is not linear to the number of CPUs used, because each thread needs to reach the correct starting line, which can takes time for very large files. (i.e. a 4.6 Gig fastq file cann contain up to 10**10 lines)

## Warning
   * The input files should be compressed using the bzip2 protocol!
   * The 4 input files should be paired: i.e. the reads are ordered similarly and match between the files

## bzip2 decompression and encoding
   * By default, this software uses the C header "bzlib.h" which provides the fastest implementation of bzip2 library. However, it is also possible to use a go bzip2 encoding/library, with the option -use_bzip2_go_lib which replaces the bzlib.h library, but is significantly slower (1.2 to 1.5 x slower).

##

## Example

```
wget http://enhancer.sdsc.edu/spreissl/Test/SP176_177_P56_I1.fastq.bz2
wget http://enhancer.sdsc.edu/spreissl/Test/SP176_177_P56_I2.fastq.bz2
wget http://enhancer.sdsc.edu/spreissl/Test/SP176_177_P56_R1.fastq.bz2
wget http://enhancer.sdsc.edu/spreissl/Test/SP176_177_P56_R1.fastq.bz2

 ATACdemultiplex -fastq_I1 SP176_177_P56_I1.fastq.bz2 -fastq_I2 SP176_177_P56_I2.fastq.bz2 -fastq_R1 SP176_177_P56_R1.fastq.bz2 -fastq_R2 SP176_177_P56_R2.fastq.bz2  --max_nb_reads 100000
```

## Performance

Depending of the options chosen, the speed is between 1.5 Meg/s (for bzip fastq files) to 2.0 Meg/s. We also prodive a multi-threading option. However, the multi-threading does lot reduce the speed linearly (each thread needs to reach their respective starting lines, which takes some time!), but rather sublinerarly.


## Contact

Olivier Poirion (PhD)
	* oporion@ucsd.edu
