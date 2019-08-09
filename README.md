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
  * Download binaries: ![https://golang.org/dl/](https://golang.org/dl/)
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