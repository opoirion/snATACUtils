/* Suite of functions dedicated to process BAM or BED files */

package main

import(
	"flag"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/biogo/hts/bgzf"
	"os"
	"log"
	"fmt"
	"bufio"
	"io"
	"strings"
	"time"
	"path"
	"sync"
	"bytes"
	"strconv"
	"sort"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
)


/*BAMFILENAME bam file name (input) */
var BAMFILENAME string

/*BEDFILENAME bam file name (input) */
var BEDFILENAME string

/*BEDFILENAMES bam file name (input) */
var BEDFILENAMES utils.ArrayFlags

/*FILENAMEOUT bam file name (output) */
var FILENAMEOUT string

/*OUTPUTDIR output directory */
var OUTPUTDIR string

/*THREADNB number of threads for reading the bam file */
var THREADNB int

/*DOWNSAMPLE number of threads for reading the bam file */
var DOWNSAMPLE float64

/*NUCLEIFILE file with cells ID */
var NUCLEIFILE string

/*NUCLEIINDEX file with cells ID */
var NUCLEIINDEX string

/*CELLIDDICT cell ID<->dict */
var CELLIDDICT map[string]bool

/*INPUTFNAMEINDEX cell ID<->dict */
var INPUTFNAMEINDEX map[string]int

/*CELLIDDICTMULTIPLE cell ID<->dict */
var CELLIDDICTMULTIPLE map[string]map[string]bool

/*OUTFILENAMELIST bam filename<->bool */
var OUTFILENAMELIST map[string]bool

/*WRITERDICT filename<->dict */
var WRITERDICT map[string]*os.File

/*BAMWRITERDICT filename<->dict */
var BAMWRITERDICT map[string]*bam.Writer

/*BEDWRITERDICT filename<->dict */
var BEDWRITERDICT map[string]io.WriteCloser

/*DIVIDE  dividing the bam file tool */
var DIVIDE bool

/*DIVIDEPARALLEL  dividing the bam file tool */
var DIVIDEPARALLEL bool

/*ADDRG  add rg group to bam file */
var ADDRG bool

/*CREATECELLINDEX  dividing the bam file tool */
var CREATECELLINDEX bool

/*BEDTOBEDGRAPHCHAN bed_to_bedgraph channel dict */
var BEDTOBEDGRAPHCHAN map[int]chan map[int]int

/*STRTOINTCHAN int -> string map chan */
var STRTOINTCHAN chan map[string]int

/*SORTFILE sort file bool */
var SORTFILE bool

/*INTCHAN int -> string map chan */
var INTCHAN chan int

/*BINSIZE bin size for bedgraph creation */
var BINSIZE int

/*BEDTOBEDGRAPH  bedgraph creation bool */
var BEDTOBEDGRAPH bool

/*DELIMITER  delimiter (string) */
var DELIMITER string

/*SPLIT split file per chromosomes */
var SPLIT bool

/*REFCHR reference chromosomes for bedgpraph */
var REFCHR map[string]int

/*REFCHRFNAME reference chromosomes for bedgpraph */
var REFCHRFNAME utils.Filename

/*NORMBEDGPRAH norm bedpgraph value */
var NORMBEDGPRAH bool


func main() {

	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `
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


`)
		 flag.PrintDefaults()
	}

	flag.StringVar(&BEDFILENAME, "bed", "", "name of the bed file")
	flag.Var(&BEDFILENAMES, "beds", "name of the bed files")
	flag.BoolVar(&BEDTOBEDGRAPH, "bed_to_bedgraph", false,
		`transform one (-bed) or multiple (use multiple -beds option) into bedgraph`)
	flag.StringVar(&BAMFILENAME, "bam", "", "name of the bam file")
	flag.Var(&REFCHRFNAME, "refchr",
		`file with reference chromosomes to use for the bedgraph creation.
 This file can also contain the maximum chromosome size, for example (chr16<tab>90668800)`)
	flag.StringVar(&FILENAMEOUT, "out", "", "name of the output file")
	flag.StringVar(&NUCLEIFILE, "cellsID", "", "file with cell IDs")
	flag.StringVar(&OUTPUTDIR, "output_dir", "", "output directory")
	flag.BoolVar(&NORMBEDGPRAH, "norm", true, `norm bedgpraph values`)
	flag.BoolVar(&CREATECELLINDEX, "create_cell_index", false, `create cell index (cell -> read Counts) for a bam or bed file`)
	flag.StringVar(&NUCLEIINDEX, "cell_index", "", "nuclei <-> output files index")
	flag.BoolVar(&DIVIDE, "divide", false, `divide the bam/bed file according to barcode file list`)
	flag.BoolVar(&DIVIDEPARALLEL, "divide_parallel", false,
		`divide the bam file according to barcode file list using a parallel version`)
	flag.BoolVar(&SORTFILE, "sort", false, "sort output file of cell index")
	flag.BoolVar(&ADDRG, "add_rg", false, "add cell ID as RG group to bam file")
	flag.BoolVar(&SPLIT, "split", false, `split file per chromosomes`)
	flag.Float64Var(&DOWNSAMPLE, "downsample", 1.0, `Downsample the number of reads from a a bed file (downsample = 1.0 is 100% and downsample = 0.0 is 0% of the reads)`)
	flag.IntVar(&THREADNB, "threads", 1, "threads concurrency for reading bam file")
	flag.IntVar(&BINSIZE, "binsize", 50, "bin size for bedgraph creation")
	flag.StringVar(&DELIMITER, "delimiter", "\t", "delimiter used")
	flag.Parse()

	if OUTPUTDIR != "" {
		FILENAMEOUT = fmt.Sprintf("%s/%s", OUTPUTDIR, FILENAMEOUT)
	}

	if BAMFILENAME == "" && BEDFILENAME == "" && len(BEDFILENAMES) == 0 {
		panic("-bam or -bed must be specified!")
	}

	tStart := time.Now()

	switch {
	case ADDRG:
		fmt.Printf("launching InsertRGTagToBamFile...\n")
		InsertRGTagToBamFile()
	case (DIVIDE || DIVIDEPARALLEL) && NUCLEIINDEX!="":
		switch {
		case BEDFILENAME != "" && DIVIDEPARALLEL:
			fmt.Printf("launching DivideMultipleBedFileParallel...\n")
			DivideMultipleBedFileParallel()
		case BEDFILENAME != "":
			fmt.Printf("launching DivideMultipleBedFile...\n")
			DivideMultipleBedFile()
		case DIVIDEPARALLEL:
			fmt.Printf("launching DivideMultipleBamFileParallel...\n")
			DivideMultipleBamFileParallel()
		default:
			fmt.Printf("launching DivideMultipleBamFile...\n")
			DivideMultipleBamFile()
		}

		for filename := range(OUTFILENAMELIST) {
			fmt.Printf("output bam file: %s written\n", filename)
		}

	case DIVIDE && BEDFILENAME != "":
		fmt.Printf("launching DivideBed...\n")
		DivideBed()
	case DIVIDE:
		fmt.Printf("launching DivideBam...\n")
		DivideBam()
	case CREATECELLINDEX:
		if FILENAMEOUT == "" {
			log.Fatal("Error -out must be provided!")
		}

		switch{
		case BAMFILENAME != "":
			fmt.Printf("launching CreateCellIndexBam...\n")
			CreateCellIndexBam()
		case BEDFILENAME != "":
			fmt.Printf("launching CreateCellIndexBed...\n")
			CreateCellIndexBed()
		default:
			log.Fatal("error! either -bam or -bed must be provided")
		}

		if SORTFILE {
			utils.SortLogfile(FILENAMEOUT, DELIMITER, "", true, false)
		}

	case BEDTOBEDGRAPH:
		if BEDFILENAME != ""{
			BEDFILENAMES = append(BEDFILENAMES, BEDFILENAME)
		}
		fmt.Printf("launching BedToBedGraphDictMultipleFile...\n")
		BedToBedGraphDictMultipleFile()
	case SPLIT:
		SplitBedPerChr()

	case DOWNSAMPLE < 1.0 && DOWNSAMPLE > 0.0:
		downSampleBedFile()
	}

	tDiff := time.Since(tStart)
	fmt.Printf("done in time: %f s \n", tDiff.Seconds())
}

func loadRefChrMap() {
	var line, chr string
	var lineSplit []string
	var chrsize int
	var err error

	if REFCHRFNAME == "" {
		return
	}

	REFCHR = make(map[string]int)

	reader, file := REFCHRFNAME.ReturnReader(0)
	defer utils.CloseFile(file)

	for reader.Scan() {
		line = reader.Text()
		line = strings.ReplaceAll(line, " ", "\t")
		lineSplit = strings.Split(line, "\t")
		chr = lineSplit[0][3:]
		REFCHR[chr] = 1

		if len(lineSplit) > 1 {
			chrsize, err = strconv.Atoi(lineSplit[1])

			if err == nil {
				REFCHR[chr] = chrsize
			}
		}
	}
}

func downSampleBedFile() {
	var line []byte
	var split [][]byte
	checkCellIndex := false
	var buffer bytes.Buffer
	var thres int

	writeToThres := false

	ext := path.Ext(BEDFILENAME)
	count := 0
	totalNbReads := 0

	if DOWNSAMPLE > 0.5 {
		thres = int(10.0 * DOWNSAMPLE)
		writeToThres = true

	} else {
		thres = int(1.0 / DOWNSAMPLE)
	}

	sep := []byte("\t")

	if FILENAMEOUT == "" {
		FILENAMEOUT = fmt.Sprintf("%s.downsample%s", BEDFILENAME[:len(BEDFILENAME) - len(ext)], ext)
	}

	if NUCLEIFILE != "" {
		loadCellIDDict(NUCLEIFILE)
		checkCellIndex = true
	}

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0)
	defer utils.CloseFile(file)

	writer := utils.ReturnWriter(FILENAMEOUT)
	defer utils.CloseFile(writer)

	for bedReader.Scan() {

		if checkCellIndex {
			line = bedReader.Bytes()
			split = bytes.Split(line, sep)
			if !CELLIDDICT[string(split[3])] {
				continue
			}
		}

		count++

		switch {

		case !writeToThres && count >= thres:
			line = bedReader.Bytes()
			buffer.Write(line)
			buffer.WriteRune('\n')
			writer.Write(buffer.Bytes())
			buffer.Reset()
			count = 0
			totalNbReads++

		case !writeToThres:

		case writeToThres && count <= thres:
			line = bedReader.Bytes()
			buffer.Write(line)
			buffer.WriteRune('\n')
			writer.Write(buffer.Bytes())
			buffer.Reset()
			totalNbReads++

		case writeToThres && count >= 10:
			count = 0
		}
	}

	fmt.Printf("Total number of reads used: %d\n", totalNbReads)
	fmt.Printf("File: %s written \n", FILENAMEOUT)
}

/*SplitBedPerChr split a bed file per chromosome */
func SplitBedPerChr() {
	var count int
	var chrID string
	var line []byte
	var split [][]byte
	var isInside bool
	var bedFname string
	var ext string
	var buffer bytes.Buffer

	checkCellIndex := false

	sep := []byte("\t")

	if FILENAMEOUT == "" {
		FILENAMEOUT = BEDFILENAME
	}

	if NUCLEIFILE != "" {
		loadCellIDDict(NUCLEIFILE)
		checkCellIndex = true
	}

	ext = path.Ext(FILENAMEOUT)

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0)

	defer utils.CloseFile(file)

	resFileDict := make(map[string]io.WriteCloser)
	resNameDict := make(map[string]string)

	for bedReader.Scan() {
		line = bedReader.Bytes()
		count++
		split = bytes.Split(line, sep)

		if checkCellIndex {
			if !CELLIDDICT[string(split[3])] {
				continue
			}
		}

		chrID = string(split[0])

		if _, isInside = resNameDict[chrID]; !isInside {
			bedFname = fmt.Sprintf("%s.%s%s",
				FILENAMEOUT[:len(FILENAMEOUT) - len(ext)], chrID, ext)
			resNameDict[chrID] = bedFname
			resFileDict[chrID] = utils.ReturnWriter(bedFname)
			defer utils.CloseFile(resFileDict[chrID])
		}

		buffer.Write(line)
		buffer.WriteRune('\n')

		resFileDict[chrID].Write(buffer.Bytes())
		buffer.Reset()
	}

	for _,f := range resNameDict {
		fmt.Printf("file %s written\n", f)
	}
}

/*CreateCellIndexBam create cell index */
func CreateCellIndexBam() {
	var record * sam.Record
	var readID string
	var count int
	var buffer bytes.Buffer

	readIndex := make(map[string]int)

	f, err := os.Open(BAMFILENAME)
	check(err)
	defer utils.CloseFile(f)
	_, err = bgzf.HasEOF(f)
	check(err)

	fOut, err := os.Create(FILENAMEOUT)

	check(err)
	defer utils.CloseFile(fOut)

	bamReader, err := bam.NewReader(f, THREADNB)
	defer utils.CloseFile(bamReader)
	check(err)

	for {
		record, err = bamReader.Read()

	errLoop:
		switch err {
		case io.EOF:
			break errLoop
		case nil:
		default:
			fmt.Printf("ERROR: %s\n",err)
			break errLoop
		}

		if record == nil {
			break
		}

		readID = strings.SplitN(record.Name, ":", 2)[0]
		readIndex[readID]++
	}


	for readID, count = range(readIndex) {
		buffer.WriteString(readID)
		buffer.WriteString(DELIMITER)
		buffer.WriteString(strconv.Itoa(count))
		buffer.WriteRune('\n')

		fOut.Write(buffer.Bytes())

		buffer.Reset()
	}
	fmt.Printf("cell index: %s created\n", FILENAMEOUT)
}


/*CreateCellIndexBed create cell index */
func CreateCellIndexBed() {
	var line string

	var readID string
	var count int
	var buffer bytes.Buffer

	readIndex := make(map[string]int)

	fOut, err := os.Create(FILENAMEOUT)

	check(err)
	defer utils.CloseFile(fOut)

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0)
	defer utils.CloseFile(file)

	for bedReader.Scan(){
		line = bedReader.Text()

		readID = strings.Split(line, "\t")[3]
		readIndex[readID]++
	}

	for readID, count = range(readIndex) {
		if readID == "" {
			continue
		}

		buffer.WriteString(readID)
		buffer.WriteString(DELIMITER)
		buffer.WriteString(strconv.Itoa(count))
		buffer.WriteRune('\n')

		fOut.Write(buffer.Bytes())

		buffer.Reset()
	}
	fmt.Printf("cell index: %s created\n", FILENAMEOUT)
}


/*InsertRGTagToBamFile insert  */
func InsertRGTagToBamFile() {
	var record * sam.Record
	var readID string
	var aux sam.Aux

	f, err := os.Open(BAMFILENAME)
	check(err)
	defer utils.CloseFile(f)
	_, err = bgzf.HasEOF(f)
	check(err)

	ext := path.Ext(BAMFILENAME)

	filenameout := fmt.Sprintf("%s.RG_added.bam", strings.TrimSuffix(BAMFILENAME, ext))

	fOut, err := os.Create(filenameout)

	check(err)
	defer utils.CloseFile(fOut)

	bamReader, err := bam.NewReader(f, THREADNB)
	check(err)
	defer utils.CloseFile(bamReader)

	header := bamReader.Header()
	bamWriter, err := bam.NewWriter(fOut, header, THREADNB)
	check(err)
	defer utils.CloseFile(bamWriter)

	loop:
	for {
		record, err = bamReader.Read()

		switch err {
		case io.EOF:
			break loop
		case nil:
		default:
			fmt.Printf("ERROR: %s\n",err)
			break loop
		}

		if record == nil {
			break loop
		}

		readID = strings.SplitN(record.Name, ":", 2)[0]

		aux = []byte(fmt.Sprintf("RGZ%s",readID))

		record.AuxFields = append(record.AuxFields, aux)

		err = bamWriter.Write(record)
		check(err)
	}

}


/*DivideMultipleBedFileParallel divide one bed file into multiple bed files in parallel */
func DivideMultipleBedFileParallel() {
	loadCellIDIndexAndBEDWriter(NUCLEIINDEX)

	for _, file := range(WRITERDICT) {
		defer utils.CloseFile(file)
	}

	for _, file := range(BEDWRITERDICT) {
		defer utils.CloseFile(file)
	}

	index:= 0

	INPUTFNAMEINDEX = make(map[string]int)

	for filename := range(OUTFILENAMELIST) {
		INPUTFNAMEINDEX[filename] = index
		index++
	}

	var waiting sync.WaitGroup
	waiting.Add(THREADNB)

	for i := 0; i < THREADNB;i++ {
		go divideMultipleBedFileOneThread(i, &waiting)
	}

	waiting.Wait()
}

/*DivideMultipleBamFileParallel divide a bam file into multiple bam files in parallel */
func DivideMultipleBamFileParallel() {
	f, err := os.Open(BAMFILENAME)
	check(err)
	defer utils.CloseFile(f)
	_, err = bgzf.HasEOF(f)
	check(err)

	bamReader, err := bam.NewReader(f, THREADNB)
	check(err)
	defer utils.CloseFile(bamReader)

	header := bamReader.Header()
	loadCellIDIndexAndBAMWriter(NUCLEIINDEX, header)
	check(err)

	for _, file := range(WRITERDICT) {
		defer utils.CloseFile(file)
	}

	for _, file := range(BAMWRITERDICT) {
		defer utils.CloseFile(file)
	}

	index:= 0

	INPUTFNAMEINDEX = make(map[string]int)

	for filename := range(OUTFILENAMELIST) {
		INPUTFNAMEINDEX[filename] = index
		index++
	}

	var waiting sync.WaitGroup
	waiting.Add(THREADNB)

	for i := 0; i < THREADNB;i++ {
		go divideMultipleBamFileOneThread(i, &waiting)
	}

	waiting.Wait()
}

func divideMultipleBamFileOneThread(threadID int, waiting *sync.WaitGroup){
	defer waiting.Done()

	var record * sam.Record
	var readID string
	var bamWriter *bam.Writer
	var isInside bool
	var filename string

	f, err := os.Open(BAMFILENAME)
	check(err)
	defer utils.CloseFile(f)
	_, err = bgzf.HasEOF(f)
	check(err)

	bamReader, err := bam.NewReader(f, THREADNB)
	check(err)
	defer utils.CloseFile(bamReader)

	count := 0
	chunkSize := len(CELLIDDICTMULTIPLE) / (THREADNB - 1)
	startIndex := chunkSize * threadID
	endIndex := chunkSize * (threadID + 1)

	loop:
	for {
		record, err = bamReader.Read()
		count++

		switch err {
		case io.EOF:
			break loop
		case nil:
		default:
			fmt.Printf("ERROR: %s\n",err)
			break loop
		}

		if record == nil {
			break loop
		}

		readID = strings.SplitN(record.Name, ":", 2)[0]

		if  _, isInside = CELLIDDICTMULTIPLE[readID];isInside {

			for filename = range(CELLIDDICTMULTIPLE[readID]) {
				if !((startIndex <= INPUTFNAMEINDEX[filename]) && (INPUTFNAMEINDEX[filename]  < endIndex)) {
					continue
				}

				bamWriter = BAMWRITERDICT[filename]
				err = bamWriter.Write(record)
				if err != nil {
					fmt.Printf("ERROR: %s\n", err)
				}
				check(err)
			}
		}
	}
}

func divideMultipleBedFileOneThread(threadID int, waiting *sync.WaitGroup){
	defer waiting.Done()

	var readID string
	var isInside bool
	var filename string
	var line string
	var buffer bytes.Buffer

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0)
	defer utils.CloseFile(file)

	count := 0
	chunkSize := len(CELLIDDICTMULTIPLE) / (THREADNB - 1)
	startIndex := chunkSize * threadID
	endIndex := chunkSize * (threadID + 1)

	for bedReader.Scan() {
		line = bedReader.Text()
		count++

		readID = strings.Split(line, "\t")[3]

		if  _, isInside = CELLIDDICTMULTIPLE[readID];isInside {
			buffer.WriteString(line)
			buffer.WriteRune('\n')

			for filename = range(CELLIDDICTMULTIPLE[readID]) {
				if !((startIndex <= INPUTFNAMEINDEX[filename]) && (INPUTFNAMEINDEX[filename]  < endIndex)) {
					continue
				}
				BEDWRITERDICT[filename].Write(buffer.Bytes())
			}
			buffer.Reset()
		}
	}
}


/*DivideMultipleBedFile divide the bam file */
func DivideMultipleBedFile() {
	var line string
	var readID string
	var isInside bool
	var filename string
	var buffer bytes.Buffer

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0)
	defer utils.CloseFile(file)

	loadCellIDIndexAndBEDWriter(NUCLEIINDEX)

	for _, file := range(WRITERDICT) {
		defer utils.CloseFile(file)
	}

	for _, file := range(BEDWRITERDICT) {
		defer utils.CloseFile(file)
	}


	count := 0

	for bedReader.Scan() {
		line = bedReader.Text()
		count++

		readID = strings.Split(line, "\t")[3]

		if  _, isInside = CELLIDDICTMULTIPLE[readID];isInside {
			buffer.WriteString(line)
			buffer.WriteRune('\n')

			for filename = range(CELLIDDICTMULTIPLE[readID]) {
				BEDWRITERDICT[filename].Write(buffer.Bytes())
			}

			buffer.Reset()
		}
	}
}

/*DivideMultipleBamFile divide the bam file */
func DivideMultipleBamFile() {
	var record * sam.Record
	var read string
	var readID string
	var bamWriter *bam.Writer
	var isInside bool
	var filename string

	f, err := os.Open(BAMFILENAME)
	check(err)
	defer utils.CloseFile(f)
	_, err = bgzf.HasEOF(f)
	check(err)

	bamReader, err := bam.NewReader(f, THREADNB)
	defer utils.CloseFile(bamReader)
	check(err)

	header := bamReader.Header()
	loadCellIDIndexAndBAMWriter(NUCLEIINDEX, header)

	for _, file := range(WRITERDICT) {
		defer utils.CloseFile(file)
	}

	for _, file := range(BAMWRITERDICT) {
		defer utils.CloseFile(file)
	}

	count := 0

	loop:
	for {
		record, err = bamReader.Read()
		count++

		switch err {
		case io.EOF:
			break loop
		case nil:
		default:
			fmt.Printf("ERROR: %s\n",err)
			break loop
		}

		if record == nil {
			break loop
		}

		read = record.String()
		readID = strings.SplitN(read, ":", 2)[0]

		if  _, isInside = CELLIDDICTMULTIPLE[readID];isInside {

			for filename = range(CELLIDDICTMULTIPLE[readID]) {
				bamWriter = BAMWRITERDICT[filename]
				err = bamWriter.Write(record)
				if err != nil {
					fmt.Printf("ERROR: %s\n", err)
				}
				check(err)
			}
		}
	}
}

/*DivideBed divide a bed file */
func DivideBed() {
	var readID string
	var line string
	var buffer bytes.Buffer

	if FILENAMEOUT == "" {
		panic("-out must be specified!")
	}

	if NUCLEIFILE == "" {
		panic("-cellsID must be specified!")
	}

	loadCellIDDict(NUCLEIFILE)

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0)
	bedWriter := utils.ReturnWriter(FILENAMEOUT)

	defer utils.CloseFile(file)
	defer utils.CloseFile(bedWriter)

	count := 0

	for bedReader.Scan() {
		line = bedReader.Text()
		count++

		readID = strings.Split(line, "\t")[3]

		if value := CELLIDDICT[readID];value {
			buffer.WriteString(line)
			buffer.WriteRune('\n')
			bedWriter.Write(buffer.Bytes())
			buffer.Reset()
		}
	}
}

/*sum ... */
func sum(array []int)int{
	sum := 0
	var val int

	for _, val = range array {
		sum += val
	}

	return sum
}


/*BedToBedGraphDictMultipleFile Transform multiple bed files given using -bed arg */
func BedToBedGraphDictMultipleFile(){
	for _, file := range BEDFILENAMES {
		if _, err := os.Stat(file); os.IsNotExist(err) {
			log.Fatal(fmt.Sprintf("### bed file: %s do not exists!!!\n", file))
		}
	}

	tStart := time.Now()
	checkCellIndex := false

	if NUCLEIFILE != "" {
		loadCellIDDict(NUCLEIFILE)
		checkCellIndex = true
	}

	INTCHAN = make(chan int, len(BEDFILENAMES))
	STRTOINTCHAN = make(chan map[string]int, len(BEDFILENAMES))
	BEDTOBEDGRAPHCHAN = make(map[int]chan map[int]int, THREADNB)

	for i := 1 ; i < 250 ; i++ {
		BEDTOBEDGRAPHCHAN[i] = make(chan map[int]int, THREADNB)
	}

	guard := make(chan struct{}, THREADNB)

	var waiting sync.WaitGroup
	waiting.Add(len(BEDFILENAMES))

	for _, file := range BEDFILENAMES {
		guard <- struct{}{}
		go bedToBedGraphDictOneThread(file, &waiting, checkCellIndex)
		<-guard
	}

	waiting.Wait()

	tDiff := time.Since(tStart)
	fmt.Printf("iterating through bed files done in: %f sec \n", tDiff.Seconds())

	collectAndProcessMultipleBedGraphDict(FILENAMEOUT)
}

func collectAndProcessMultipleBedGraphDict(filenameout string) {
	listNbReads := extractIntfromChan(INTCHAN)
	chroHashDict := extractStrIntDictFromChan(STRTOINTCHAN)
	totNbReads := sum(listNbReads)
	guard := make(chan struct{}, THREADNB)

	fmt.Printf("sumf of number READS: %d\n", totNbReads)

	var waitingSort sync.WaitGroup
	var scale float64

	loadRefChrMap()

	if filenameout == "" {
		filenameout = BEDFILENAMES[0]
	}

	ext := path.Ext(filenameout)

	if ext == ".bedgraph" {
		filenameout = filenameout[0:len(filenameout) - len(ext)]
	}

	if NORMBEDGPRAH {
		scale = (float64(totNbReads) / 1e6) * (float64(BINSIZE) / 1e3)

	} else {
		scale = 1.0
	}

	tStart := time.Now()

	chroList := []string{}
	fileList := []string{}

	for chro := range chroHashDict {
		if len(chroHashDict) > 0 {
			chroList = append(chroList, chro)
		}
	}
	sort.Strings(chroList)

	filterChro := len(REFCHR) > 0

	for _, chro := range chroList {
		if filterChro && REFCHR[chro] == 0 {
			continue
		}

		guard <- struct{}{}
		chroID := chroHashDict[chro]
		bedFname := fmt.Sprintf("%s.chr%s.bedgraph", filenameout, chro)
		fileList = append(fileList, bedFname)
		waitingSort.Add(1)
		go writeIndividualChrBedGraph(
			bedFname, chroID, chro, scale, REFCHR[chro], &waitingSort)
		<-guard
	}

	waitingSort.Wait()
	tDiff := time.Since(tStart)
	fmt.Printf(" writing done in: %f sec \n", tDiff.Seconds())

	if !SPLIT {
		tStart = time.Now()
		cmd := fmt.Sprintf("cat %s > %s.bedgraph", strings.Join(fileList," "), filenameout)

		utils.ExceCmd(cmd)
		fmt.Printf("%s.bedgraph created!\n", filenameout)

		cmd = fmt.Sprintf("rm %s", strings.Join(fileList, " "))
		utils.ExceCmd(cmd)
		fmt.Printf(" concatenation and cleaning done in: %f sec \n", tDiff.Seconds())
	} else {
		for _,f := range fileList {
			fmt.Printf("file: %s created\n", f)
		}
	}
}


/*writeIndividualChrBedGraph write an  bedgraph file for a unique chromosome */
func writeIndividualChrBedGraph(bedFname string, chroID int, chro string,
	scale float64, limit int, waiting * sync.WaitGroup){
	defer waiting.Done()

	fWrite, err := os.Create(bedFname)
	check(err)
	defer utils.CloseFile(fWrite)

	bedgraphDict := extractIntDictFromChan(BEDTOBEDGRAPHCHAN[chroID])
	var buffer bytes.Buffer

	chroPosList := make([]int, 0, len(bedgraphDict))

	for key := range bedgraphDict {
		chroPosList = append(chroPosList, key)
	}

	sort.Ints(chroPosList)
	var pos, posEnd int
	var value float64
	currentPos := 0

	buffSize := 0

	for _, binnb := range chroPosList {
		pos = binnb * BINSIZE

		if pos != currentPos {
			buffer.WriteString("chr")
			buffer.WriteString(chro)
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.Itoa(currentPos))
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.Itoa(pos))
			buffer.WriteRune('\t')
			buffer.WriteRune('0')
			buffer.WriteRune('\n')

			buffSize++

			if buffSize > 50000 {
				_, err = fWrite.Write(buffer.Bytes())
				check(err)
				buffer.Reset()
				buffSize = 0

			}
		}

		value = float64(bedgraphDict[binnb]) / scale

		posEnd = pos + BINSIZE

		if limit > 1 {
			if   pos > limit {
				fmt.Printf(
					"!!!! Warning starting genomic position: %d greater than chr%s limit: %d. Skipping...\n",
					pos, chro, limit)
				break

			}
			if   posEnd > limit {
				posEnd = limit
			}
		}

		buffer.WriteString("chr")
		buffer.WriteString(chro)
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.Itoa(pos))
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.Itoa(posEnd))
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.FormatFloat(value, 'f', 2, 32))
		buffer.WriteRune('\n')

		buffSize++

		if buffSize > 50000 {
			_, err = fWrite.Write(buffer.Bytes())
			check(err)
			buffer.Reset()
			buffSize = 0

		}

		currentPos = pos + BINSIZE
	}

	_, err = fWrite.Write(buffer.Bytes())
	check(err)
	buffer.Reset()
}

/*extractIntfromChan Transform multiple bed files given using -bed arg */
func extractIntfromChan(channel chan int) ([]int) {

	list := []int{}

	loop:
	for {
		select{
		case val := <-channel:
			list = append(list, val)
		default:
			break loop
		}
	}

	return list
}

/*extractStrIntDictFromChan extract and update multiple string -> int dict */
func extractStrIntDictFromChan(channel chan map[string]int) (map[string]int) {

	dict := make(map[string]int)
	var isInside bool

	loop:
	for {
		select{
		case statsDict := <-channel:
			for key, value := range statsDict {
				if _, isInside = dict[key];!isInside {
					dict[key] = value
				}
			}
		default:
			break loop
		}
	}

	return dict
}


/*extractIntDictFromChan Transform multiple bed files given using -bed arg */
func extractIntDictFromChan(channel chan map[int]int) (map[int]int) {

	dict := make(map[int]int)

	loop:
	for {
		select{
		case statsDict := <-channel:
			for key, value := range statsDict {
				dict[key] += value
			}
		default:
			break loop
		}
	}

	return dict
}

/*bedToBedGraphDictOneThread Transform a bed file into a bedgraph dict */
func bedToBedGraphDictOneThread(bed string, waiting *sync.WaitGroup, checkCellIndex bool){
	defer waiting.Done()
   	bedReader, file := utils.ReturnReader(bed, 0)
	defer utils.CloseFile(file)

	var split []string
	var chroStr string
	var chroIndex int
	var pos int
	var lineNb int
	var err error
	var nbReads int
	var line string
	var isInside bool

	bedtobedgraphdict := make(map[int]map[int]int)
	chrDict := make(map[string]int)

	for i := 1; i < 250; i++ {
		bedtobedgraphdict[i] = make(map[int]int)
	}

	uncommon := 200

	for bedReader.Scan() {
		line = bedReader.Text()

		split = strings.Split(line, "\t")

		if checkCellIndex {
			if isInside = CELLIDDICT[split[3]];!isInside {
				continue
			}
		}

		nbReads++

		chroStr = strings.TrimPrefix(split[0], "chr")
		chroIndex, err = strconv.Atoi(chroStr)

		if err != nil {

			if _, isInside = chrDict[chroStr];!isInside {

				if len(chroStr) == 1 {
					chrDict[chroStr] = int(chroStr[0]) + 100
				} else {
					chrDict[chroStr] = uncommon
					uncommon++

					if _, isInside = bedtobedgraphdict[uncommon];!isInside {
						bedtobedgraphdict[uncommon] = make(map[int]int)

					}
				}
			}

			chroIndex = chrDict[chroStr]
		}

		pos, err = strconv.Atoi(split[1])

		if err != nil {
			log.Fatal(fmt.Sprintf("Error with line %s at position: %d\n",
				line, lineNb))
		}

		bedtobedgraphdict[chroIndex][pos / BINSIZE]++
	}

	INTCHAN <- nbReads
	STRTOINTCHAN <- chrDict

	for key := range bedtobedgraphdict {



		if len(bedtobedgraphdict[key]) > 0 {
			if _, isInside = BEDTOBEDGRAPHCHAN[key];!isInside {
				BEDTOBEDGRAPHCHAN[key] = make(chan map[int]int, THREADNB)
			}

			if key < 50 {
				chrDict[strconv.Itoa(key)] = key
			}

			BEDTOBEDGRAPHCHAN[key] <- bedtobedgraphdict[key]
		}
	}
}


/*DivideBam divide the bam file */
func DivideBam() {
	var record * sam.Record
	var read string
	var readID string

	if FILENAMEOUT == "" {
		panic("-out must be specified!")
	}

	if NUCLEIFILE == "" {
		panic("-cellsID must be specified!")
	}

	loadCellIDDict(NUCLEIFILE)

	f, err := os.Open(BAMFILENAME)
	check(err)
	defer utils.CloseFile(f)
	_, err = bgzf.HasEOF(f)
	check(err)

	fWrite, err := os.Create(FILENAMEOUT)
	check(err)
	defer utils.CloseFile(fWrite)

	bamReader, err := bam.NewReader(f, THREADNB)
	check(err)
	defer utils.CloseFile(bamReader)

	header := bamReader.Header()
	bamWriter, err := bam.NewWriter(fWrite, header, THREADNB)
	check(err)
	defer utils.CloseFile(bamWriter)
	count := 0

	loop:
	for {
		record, err = bamReader.Read()
		count++

		switch err {
		case io.EOF:
			break loop
		case nil:
		default:
			fmt.Printf("ERROR: %s\n",err)
			break loop
		}

		if record == nil {
			break loop
		}

		read = record.String()
		readID = strings.SplitN(read, ":", 2)[0]

		if value := CELLIDDICT[readID];value {
			err = bamWriter.Write(record)
			if err != nil {
				fmt.Printf("ERROR: %s\n", err)
			}
			check(err)
		}
	}
}


func loadCellIDDict(fname string) {
	f, err := os.Open(fname)
	check(err)
	defer utils.CloseFile(f)
	scanner := bufio.NewScanner(f)

	CELLIDDICT = make(map[string]bool)

	for scanner.Scan() {
		line := scanner.Text()
		line = strings.Split(line, "\t")[0]
		line = strings.Trim(line, " \t")
		CELLIDDICT[line] = true
	}
}

func loadCellIDIndexAndBAMWriter(fname string, header *sam.Header) {
	f, err := os.Open(fname)
	check(err)
	defer utils.CloseFile(f)
	scanner := bufio.NewScanner(f)
	var filePath string

	CELLIDDICTMULTIPLE = make(map[string]map[string]bool)
	WRITERDICT = make(map[string]*os.File)
	BAMWRITERDICT = make(map[string]*bam.Writer)
	OUTFILENAMELIST = make(map[string]bool)

	for scanner.Scan() {
		line := scanner.Text()
		split := strings.Split(line, "\t")
		cellid := split[0]
		filename := split[1]

		ext := path.Ext(filename)

		if ext != ".bam" {
			filePath = fmt.Sprintf("%s/%s.bam", OUTPUTDIR, filename)
		} else {
			filePath = fmt.Sprintf("%s/%s", OUTPUTDIR, filename)
		}

		if _, isInside := WRITERDICT[filename]; !isInside {
			WRITERDICT[filename], err = os.Create(filePath)
			check(err)

			BAMWRITERDICT[filename], err = bam.NewWriter(WRITERDICT[filename], header, THREADNB)
			check(err)
			OUTFILENAMELIST[filename] = true
		}

		if _, isInside := CELLIDDICTMULTIPLE[cellid]; !isInside {
			CELLIDDICTMULTIPLE[cellid] = make(map[string]bool)
		}

		CELLIDDICTMULTIPLE[cellid][filename] = true
	}
}


func loadCellIDIndexAndBEDWriter(fname string) {
	var filePath string
	f, err := os.Open(fname)
	check(err)
	defer utils.CloseFile(f)
	scanner := bufio.NewScanner(f)

	CELLIDDICTMULTIPLE = make(map[string]map[string]bool)
	BEDWRITERDICT = make(map[string]io.WriteCloser)
	OUTFILENAMELIST = make(map[string]bool)

	for scanner.Scan() {
		line := scanner.Text()
		split := strings.Split(line, "\t")
		cellid := split[0]
		filename := split[1]

		ext := path.Ext(filename)

		if ext != ".gz" && ext != ".bed"  {
			filePath = fmt.Sprintf("%s/%s.bed.gz", OUTPUTDIR, filename)
		} else {
			filePath = fmt.Sprintf("%s/%s", OUTPUTDIR, filename)
		}

		if _, isInside := BEDWRITERDICT[filename]; !isInside {
			BEDWRITERDICT[filename] = utils.ReturnWriter(filePath)
			check(err)
			OUTFILENAMELIST[filename] = true
		}

		if _, isInside := CELLIDDICTMULTIPLE[cellid]; !isInside {
			CELLIDDICTMULTIPLE[cellid] = make(map[string]bool)
		}

		CELLIDDICTMULTIPLE[cellid][filename] = true
	}
}

func check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
