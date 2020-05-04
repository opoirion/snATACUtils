# single-cell ATAC-seq tools

Suite of tools to maniputate files related to single-cell ATAC experiments, such as fastq files or single cell bed files (below an example)

```bash
# Example of BED file related to single-cell ATAC-Seq. The last column designates the cell barcode
hr17   14066243        14066440        AACGAGAGCTAAACCCGAGATA
chr19   42120584        42120781        AACGAGAGCTAAACCCGAGATA
chr19   42120603        42120799        AACGAGAGCTAAACCCGAGATA
chr6    36835546        36835742        AACGAGAGCTAAACCCGAGATA
chr17   14066259        14066455        AACGAGAGCTAAACCCGAGATA
chr16   79321071        79321267        AACGAGAGCTAAACCCGAGATA
chr19   45396965        45397161        AACGAGAGCTAAACCCGAGATA
```

## Installation

* Install and configure a golang compiler (if not existing)
  * Download binaries: [https://golang.org/dl/](https://golang.org/dl/)
  * Configure $GOPATH/$GOBIN

```bash
	#In .bashrc or .zshrc
	export GOROOT=$HOME/go # or wherever is you go folder
	export GOBIN=$HOME/go/local/bin # or wherever is your local bin folder for go exectuable
	export GOPATH=$HOME/go/code/:$HOME/code

	PATH=$GOPATH:$GOROOT:$PATH
	PATH=$HOME/go/bin/:$GOBIN:$PATH

```
  * source your init file `source ~/.bashrc`

* Install the package

```bash
	go get -v -u gitlab.com/Grouumf/ATACdemultiplex/...
```

The ATAC tools executable are located in your `$GOBIN` folder and should be in your global path

```bash
ATACdemultiplex -h
ATACCellTSS -h
ATACeQTLUtils -h
ATACMatUtils -h
ATACSimUtils -h
ATACtools -h
ATACTopFeatures -h
BAMutils -h
ATACAnnotateRegions -h
```

## ATACdemultiplex: Fastq files demultiplexification
Tools to insert snATAC-Seq barcodes from multiple files inside read IDs to create new fastq files containing the cell ID barcode at the begining of each read. See [https://gitlab.com/Grouumf/ATACdemultiplex/tree/master/ATACdemultiplex](https://gitlab.com/Grouumf/ATACdemultiplex/tree/master/ATACdemultiplex)

## ATACCellTSS: Computing cell / cluster TSS

```bash
USAGE: ATACCellTSS -bed <filename> -ygi/tss <filename> -xgi <filename>
                    (optionnal -out <string> -flank <int> -smoothing <int> -boundary <int> -cluster <filename> -flank_size).

if -cluster is provided, TSS is computed per cluster and -xgi argument is ignored. THe cluster file should contain cluster and cell ID with the following structure for each line: clusterID<TAB>cellID
```

## ATACMatUtils: Suite of functions dedicated to analyze intersection with genomic regions from a peak file (in bed format) to create a sparse matrix (cell x genomic regions)

```bash
#################### MODULE TO CREATE (cell x genomic region) SPARSE MATRIX ########################
"""Boolean / interger Peak matrix """
transform one (-bed) or multiple bed files into a sparse matrix
USAGE: ATACMatTools -coo -bed  <bedFile> -ygi <bedFile> -xgi <fname> (-threads <int> -out <fname> -use_count -taiji -bed <fileName2>)

"""Create a cell x bin matrix: -bin """
transform one (-bed) or multiple (use multiple -bed options) bed file into a bin (using float) sparse matrix. If ygi provided, reads intersecting these bin are ignored

USAGE: ATACMatTools -bin -bed  <bedFile> (optional -ygi <bedFile> -xgi <fname> -bin_size <int> -ygi_out <string> -norm -taiji -coo)

"""Count the number of reads in peaks for each cell: -count """
USAGE: ATACMatTools -count  -xgi <fname> -ygi <bedfile> -bed <bedFile> (optionnal: -out <fname> -norm)

"""Merge multiple matrices results into one output file: -merge """
It can be used to convert taiji to coo or coo to taiji formats.
USAGE: ATACMatTools -coo/taiji -merge -xgi <fname> -in <matrixFile1> -in <matrixFile2> ... (optional -bin -use_count -out <fname>)
```

## BAMutils: Suite of functions dedicated to process BAM or BED files

```bash
#################### Suite of functions dedicated to process BAM or BED files ########################

-bed_to_bedgraph: Transform one (-bed) or multiple (use multiple -beds option) into bedgraph
USAGE: BAMutils -bed_to_bedgraph -bed <fname> (-out <fname> -threads <int> -cellsID <fname> -split -binsize <int> -refchr <filename>)

-create_cell_index: Create cell index (cell -> read Counts) for a bam or bed file
USAGE: BAMutils -create_cell_index -bed/bam <name> -out <output name> (-sort)

-divide: Divide the bam/bed file according to barcode file list
USAGE: BAMutils -divide -bed/bam <fname> (-cell_index <fname> -threads <int> -cellsID <fname> -out <fname>)

-divide_parallel: Divide the bam file according to barcode file list using a parallel version
USAGE: BAMutils -divide_parallel -cell_index <fname> -bed/bam <fname> (-threads <int>)

-split: Split file per chromosomes
USAGE: BAMutils -split -bed <bedfile> (-out <string> -cellsID <string>)

-downsample: Downsample the number of reads from a a bed file (downsample = 1.0 is 100 perc. and downsample = 0.0 is 0 perc. of the reads)
USAGE: BAMutils -downsample <float> -bed <bedfile> (-out <string> -cellsID <string>)

-bamtobed: Transform a 10x BAM file to a bed file with each read in a new line and using the "CB:Z" field as barcode
USAGE: BAMutils -bamtobed -bam <filename> -out <bedfile> (-optionnal -cellsID <filename> -threads <int> -tag <string>)
```

## ATACTopFeatures: Module to inter significant cluster peaks using a peak list, a bed file and cell ID <-> cluster ID file


```bash
#################### MODULE TO INFER SIGNIFICANT CLUSTER PEAKS ########################

"""full individual chi2 computation for each peak with FDR correction using Benjamini-Hochberg correction. Not recommended because using golang suboptimal chi2 implementation"""
USAGE: ATACTopFeatures -chi2 -bed <fname> -peak <fname> -cluster <fname> (optionnal -out <string> -threads <int> -alpha <float> -write_all -split <int>)

"""Create contingency table for each feature and each cluster"""
USAGE: ATACTopFeatures -create_contingency -bed <fname> -peak <fname> -cluster <fname> (optionnal -out <string> -threads <int>)

"""correct feature pvalue for multiple tests performed or each cluster"""
USAGE: ATACTopFeatures -pvalue_correction -ptable <fname> (optionnal -out <string> -threads <int> -alpha <float> -write_all)
```

* Once the contingency table is created, it is preferable to use Python (or R) to infer the p-values. We wrote a python script to handle the contingency table using multithreading and the scipy package here: [https://gitlab.com/Grouumf/ATACdemultiplex/blob/master/scripts/snATAC_feature_selection](https://gitlab.com/Grouumf/ATACdemultiplex/blob/master/scripts/snATAC_feature_selection)


## ATACSimUtils: Suite of functions dedicated to generate Simulated snATAC-Seq data

```bash
#################### MODULE TO CREATE SIMULATED SINGLE CELL ATAC BED FILE ########################

USAGE: ATACSimUtils -simulate -nb <int> -mean <float> std <float> -bed <bedfile> (-threads <int> -out <string> -tag <string>)
```

## ATACtools: Suite of functions dedicated to pre/post process generic files related to snATAC pipeline

```bash
#################### Suite of functions dedicated to pre/post process files related to snATAC pipeline ########################

-bed_to_cicero: format a bed to cicero input (ex: chr1\t1215523\t1216200\tcellID -> chr1_1215523_121620,tcellID,1)
USAGE: ATACtools -bed_to_cicero -filename <bedfile> (-filenames <bedfile2> -filenames <bedfile3> ...  -ignoreerror)

-create_ref_fastq: Create a ref FASTQ file using a reference barcode list
USAGE: ATACtools -create_ref_bed -filename <fname> (-ref_barcode_list <fname> -tag <string> -output <string>)

-merge: Merge input log files together
USAGE: ATACtools -merge -filenames <fname1> -filenames <fname2> -filenames <fname3> ...  (-sortfile -delimiter "<string>" -ignoreerror -ignore_sorting_category)

-sortfile: Sort key -> values (i.e.: <key><SEP><value>) file
USAGE: ATACtools -sortfile -filename <fname> (-delimiter <string> -ignoreerror -ignore_sorting_category)

-write_compl: Write the barcode complement of a fastq files
USAGE: ATACtools -write_compl <fastq_file> (-compl_strategy <"split_10_compl_second"/"split_10_compl_first"> -tag <string>)

-scan: Scan a file and determine the number of line
USAGE ATACtools -scan -filename <string> (-printlastline -printlastlines <int> -search_in_line <string> -gotoline <int>)

-create_barcode_dict: Create a barcode key / value count file
USAGE: ATACtools -create_barcode_dict -filename <fname> (-sortfile -delimiter <string>)

-clean: clean file from unwanted lines
USAGE: ATACtools -clean -filename <fname> -output filename -clean_pattern "\n"
```


## ATACAnnotateregions: Module to annotate genomic regions from bed file using a reference bed file containing annotation

```bash

#################### MODULE TO ANNOTATE GENOMIC REGIONS FROM BED FILES ########################
This software presents some similarities with bedtools usage however it provides better customisations for bed file annotation when comparing two bed files with interesecting regions

"""Annotate bed file using a reference bed file containing the annotations"""
USAGE: ATACAnnotateRegions -bed <file> -ref <file> (optionnal -out <string> -unique -unique_ref -intersect -write_ref -edit -ref_sep "[3]int" -ref_symbol "[]int|str" -diff -stdout -annotate_line)

for -ref_sep and -ref_symbol options, the input should be a string of numbers separated by whitespace and delimited with ". -ref_sep needs exactly three positions: 1) for the chromosomes column, 2) for the begining and 3) for the end of the region

Example: ATACAnnotateRegions -bed regionToAnnotate.bed -ref referenceAnnotation.tsv -ref_sep "0 1 2" -ref_symbol "4 5"

-ref_sep and -symbol_pos can be blank with " or comma separated: i.e."0 1 2" or 0,1,2. Also symbol_pos an be a generic string to annotate all ref regions

Here the three first columns of referenceAnnotation.tsv will be used to identify chromosome (column 0), start (column 1), and end (column 2) of each region, and regionToAnnotate.bed will be annotatd using columns 4 and 5 from referenceAnnotation.tsv

```

## ATACeQTLUtils: Module to deal with eQTL from bed files and snATAC-Seq

In construction. See `ATACeQTLUtils -h`
