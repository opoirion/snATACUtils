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
ATACdemultiplex -hUsage of ATACdemultiplex:
  -fastq_I1 string
        fastq index file index paired read 1
  -fastq_I2 string
        fastq index file index paired read 2
  -fastq_R1 string
        fastq read file index paired read 1
  -fastq_R2 string
        fastq read file index paired read 2
  -max_nb_reads int
        <OPTIONAL> max number of reads to process (default 0 => None)
  -nbThreads int
        number of threads to use (default 1)
  -taglength int
        <OPTIONAL> number of nucleotides to consider at the end
         and begining (default 8) (default 8)

```

## multithreads

The speed is not linear to the number of CPUs used, because each thread needs to reach the correct starting line, which can takes time for very large files. (i.e. a 4.6 Gig fastq file cann contain up to 10**10 lines)

## Warning
   * The input files should be compressed using the bzip2 protocol!
   * The 4 input files should be paired: i.e. the reads are ordered similarly and match between the files

## Example

```
wget http://enhancer.sdsc.edu/spreissl/Test/SP176_177_P56_I1.fastq.bz2
wget http://enhancer.sdsc.edu/spreissl/Test/SP176_177_P56_I2.fastq.bz2
wget http://enhancer.sdsc.edu/spreissl/Test/SP176_177_P56_R1.fastq.bz2
wget http://enhancer.sdsc.edu/spreissl/Test/SP176_177_P56_R1.fastq.bz2

 ATACdemultiplex -fastq_I1 SP176_177_P56_I1.fastq.bz2 -fastq_I2 SP176_177_P56_I2.fastq.bz2 -fastq_R1 SP176_177_P56_R1.fastq.bz2 -fastq_R2 SP176_177_P56_R2.fastq.bz2  --max_nb_reads 100000
```

## Performance

We have currently a speed of 10s for 100K reads


## Contact

Olivier Poirion (PhD)
	* oporion@ucsd.edu
