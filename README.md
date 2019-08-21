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
```

##ATACdemultiplex: Fastq files demultiplexification
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
"""Boolean peak matrix: -coo """
transform one (-bed) or multiple (use multiple -beds option) into a boolean sparse matrix in COO format
                USAGE: ATACMatTools -coo -bed  <bedFile> -ygi <bedFile> -xgi <fname>

"""Create a cell x bin matrix: -bin """
transform one (-bed) or multiple (use multiple -beds option) into a bin (using float) sparse matrix in COO format. If ygi provided, reads intersecting these bin are ignored

USAGE: ATACMatTools -bin -bed  <bedFile> (optionnal -ygi <bedFile> -xgi <fname>) -norm

"""Count the number of reads in peaks for each cell: -count """
USAGE: ATACMatTools -count  -xgi <fname> -ygi <bedfile> -bed <bedFile>

"""Merge multiple matrices results into one output file: -merge """
USAGE: ATACMatTools -coo -merge -xgi <fname> -in <matrixFile1> -in <matrixFile2> ...
```

## BAMutils: Suite of functions dedicated to process BAM or BED files

```bash
#################### Suite of functions dedicated to process BAM or BED files ########################

-bed_to_bedgraph: Transform one (-bed) or multiple (use multiple -beds option) into bedgraph
USAGE: BAMutils -bed_to_bedgraph -bed <fname> (-out <fname> -threads <int> -cellsID <fname> -split)

-create_cell_index: Create cell index (cell -> read Counts) for a bam or bed file
USAGE: BAMutils -create_cell_index -bed/bam <name> -out <output name> (-sort)

-divide: Divide the bam/bed file according to barcode file list
USAGE: BAMutils -divide -bed/bam <fname> (-cell_index <fname> -threads <int> -cellsID <fname>)

-divide_parallel: Divide the bam file according to barcode file list using a parallel version
USAGE: BAMutils -divide_parallel -cell_index <fname> -bed/bam <fname> (-threads <int>)

-split: Split file per chromosomes
USAGE: BAMutils -split -bed <bedfile> (-out <string> -cellsID <string>)

-downsample: Downsample the number of reads from a a bed file (downsample = 1.0 is 100 perc. and downsample = 0.0 is 0 perc. of the reads)
USAGE: BAMutils -downsample <float> -bed <bedfile> (-out <string> -cellsID <string>)
```

## ATACTopFeatures: Module to inter significant cluster peaks using a peak list, a bed file and cell ID <-> cluster ID file


```bash
#################### MODULE TO INFER SIGNIFICANT CLUSTER PEAKS ########################

-chi2: Full individual chi2 computation for each peak with FDR correction using Benjamini-Hochberg correction. Not recommended because using golang suboptimal chi2 implementation
USAGE: ATACTopFeatures -chi2 -bed <fname> -peak <fname> -cluster <fname> (optionnal -out <string> -threads <int> -alpha <float> -write_all -split <int>)

-create_contingency: Create contingency table for each feature and each cluster
USAGE: ATACTopFeatures -create_contingency -bed <fname> -peak <fname> -cluster <fname> (optionnal -out <string> -threads <int>)

-pvalue_correction: Correct feature pvalue for multiple tests performed or each cluster using the Benjamini-Hochberg correction
USAGE: ATACTopFeatures -pvalue_correction -ptable <fname> (optionnal -out <string> -threads <int> -alpha <float> -write_all)
```

* Once the contingency table is created, it is preferable to use Python (or R) to infer the p-values. We wrote a python script to handle the contingency table using multithreading and the scipy package here: [https://gitlab.com/Grouumf/ATACdemultiplex/blob/master/scripts/snATAC_feature_selection] (https://gitlab.com/Grouumf/ATACdemultiplex/blob/master/scripts/snATAC_feature_selection)


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
USAGE: ATACtools -create_ref_bed -filename <fname> (-ref_barcode_list <fname> -tag <string>)

-merge: Merge input log files together
USAGE: ATACtools -merge -filenames <fname1> -filenames <fname2> -filenames <fname3> ...  (-sortfile -delimiter "<string>" -ignoreerror -ignore_sorting_category)

-sortfile: Sort key -> values (i.e.: <key><SEP><value>) file
USAGE: ATACtools -sortfile -filename <fname> (-delimiter <string> -ignoreerror -ignore_sorting_category)

-write_compl: Write the barcode complement of a fastq files
USAGE: ATACtools -write_compl <fastq_file> (-compl_strategy <"split_10_compl_second"/"split_10_compl_first"> -tag <string>)

-scan: Scan a file and determine the number of line
USAGE ATACtools -scan -filename <string> (-printlastline -printlastlines <int> -search_in_line <string> -gotoline <int>)

-create_barcode_dict: Create a barcode key / value count file
USAGE: ATACtools -create_barcode_list -filename <fname> (-sortfile -delimiter <string>)
```

## ATACeQTLUtils: Module to deal with eQTL from bed files and snATAC-Seq

In construction. See `ATACeQTLUtils -h`