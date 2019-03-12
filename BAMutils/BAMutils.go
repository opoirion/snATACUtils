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
	utils "ATACdemultiplex/ATACdemultiplexUtils"
)


/*BAMFILENAME bam file name (input) */
var BAMFILENAME string

/*BEDFILENAME bam file name (input) */
var BEDFILENAME string

/*FILENAMEOUT bam file name (output) */
var FILENAMEOUT string

/*OUTPUTDIR output directory */
var OUTPUTDIR string

/*THREADNB number of threads for reading the bam file */
var THREADNB int

/*NUCLEIFILE file with cells ID */
var NUCLEIFILE string

/*NUCLEIINDEX file with cells ID */
var NUCLEIINDEX string

/*CELLIDDICT cell ID<->dict */
var CELLIDDICT map[string]bool

/*CELLIDDINDEX cell ID<->dict */
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

func main() {
	flag.StringVar(&BEDFILENAME, "bed", "", "name of the bed file")
	flag.StringVar(&BAMFILENAME, "bam", "", "name of the bam file")
	flag.StringVar(&FILENAMEOUT, "out", "", "name of the output file")
	flag.StringVar(&NUCLEIFILE, "cellsID", "", "file with cell IDs")
	flag.StringVar(&OUTPUTDIR, "output_dir", "", "output directory")
	flag.BoolVar(&CREATECELLINDEX, "create_cell_index", false, "create cell index")
	flag.StringVar(&NUCLEIINDEX, "cell_index", "", "nuclei <-> output files index")
	flag.BoolVar(&DIVIDE, "divide", false, "divide the bam file according to barcode file list")
	flag.BoolVar(&DIVIDEPARALLEL, "divide_parallel", false, "divide the bam file according to barcode file list usng a parallel version")
	flag.BoolVar(&ADDRG, "add_rg", false, "add cell ID as RG group to bam file")
	flag.IntVar(&THREADNB, "threads", 1, "threads concurrency for reading bam file")
	flag.Parse()

	if OUTPUTDIR != "" {
		FILENAMEOUT = fmt.Sprintf("%s/%s", OUTPUTDIR, FILENAMEOUT)
	}

	if BAMFILENAME == "" && BEDFILENAME == "" {
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

		for filename, _ := range(OUTFILENAMELIST) {
			fmt.Printf("output bam file: %s written\n", filename)
		}

	case DIVIDE && BEDFILENAME != "":
		fmt.Printf("launching DivideBed...\n")
		DivideBed()
	case DIVIDE:
		fmt.Printf("launching DivideBam...\n")
		DivideBam()
	case CREATECELLINDEX:
		fmt.Printf("launching CreateCellIndex...\n")
		CreateCellIndex()
	}

	tDiff := time.Now().Sub(tStart)
	fmt.Printf("done in time: %f s \n", tDiff.Seconds())
}

/*CreateCellIndex create cell index */
func CreateCellIndex() {
	var record * sam.Record
	var readID string
	var count int
	var buffer bytes.Buffer

	readIndex := make(map[string]int)

	f, err := os.Open(BAMFILENAME)
	check(err)
	defer f.Close()
	_, err = bgzf.HasEOF(f)
	check(err)

	fOut, err := os.Create(FILENAMEOUT)

	check(err)
	defer fOut.Close()

	bamReader, err := bam.NewReader(f, THREADNB)
	defer bamReader.Close()
	check(err)

	for {
		record, err = bamReader.Read()

		switch err {
		case io.EOF:
			break
		case nil:
		default:
			fmt.Printf("ERROR: %s\n",err)
			break
		}

		if record == nil {
			break
		}

		readID = strings.SplitN(record.Name, ":", 2)[0]
		readIndex[readID]++
	}


	for readID, count = range(readIndex) {
		buffer.WriteString(readID)
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.Itoa(count))
		buffer.WriteRune('\n')

		fOut.Write(buffer.Bytes())

		buffer.Reset()
	}
}


/*InsertRGTagToBamFile insert  */
func InsertRGTagToBamFile() {
	var record * sam.Record
	var readID string
	var aux sam.Aux

	f, err := os.Open(BAMFILENAME)
	check(err)
	defer f.Close()
	_, err = bgzf.HasEOF(f)
	check(err)

	ext := path.Ext(BAMFILENAME)

	filenameout := fmt.Sprintf("%s.RG_added.bam", strings.TrimSuffix(BAMFILENAME, ext))

	fOut, err := os.Create(filenameout)

	check(err)
	defer fOut.Close()

	bamReader, err := bam.NewReader(f, THREADNB)
	defer bamReader.Close()
	check(err)

	header := bamReader.Header()
	bamWriter, err := bam.NewWriter(fOut, header, THREADNB)
	check(err)
	defer bamWriter.Close()

	for {
		record, err = bamReader.Read()

		switch err {
		case io.EOF:
			break
		case nil:
		default:
			fmt.Printf("ERROR: %s\n",err)
			break
		}

		if record == nil {
			break
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
		defer file.Close()
	}

	for _, file := range(BEDWRITERDICT) {
		defer file.Close()
	}

	index:= 0

	INPUTFNAMEINDEX = make(map[string]int)

	for filename, _ := range(OUTFILENAMELIST) {
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
	defer f.Close()
	_, err = bgzf.HasEOF(f)
	check(err)

	bamReader, err := bam.NewReader(f, THREADNB)
	defer bamReader.Close()
	check(err)

	header := bamReader.Header()
	loadCellIDIndexAndBAMWriter(NUCLEIINDEX, header)
	check(err)

	for _, file := range(WRITERDICT) {
		defer file.Close()
	}

	for _, file := range(BAMWRITERDICT) {
		defer file.Close()
	}

	index:= 0

	INPUTFNAMEINDEX = make(map[string]int)

	for filename, _ := range(OUTFILENAMELIST) {
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
	defer f.Close()
	_, err = bgzf.HasEOF(f)
	check(err)

	bamReader, err := bam.NewReader(f, THREADNB)
	defer bamReader.Close()
	check(err)

	count := 0
	chunkSize := len(CELLIDDICTMULTIPLE) / (THREADNB - 1)
	startIndex := chunkSize * threadID
	endIndex := chunkSize * (threadID + 1)

	for {
		record, err = bamReader.Read()
		count++

		switch err {
		case io.EOF:
			break
		case nil:
		default:
			fmt.Printf("ERROR: %s\n",err)
			break
		}

		if record == nil {
			break
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

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0, false)
	defer file.Close()

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

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0, false)
	defer file.Close()

	loadCellIDIndexAndBEDWriter(NUCLEIINDEX)

	for _, file := range(WRITERDICT) {
		defer file.Close()
	}

	for _, file := range(BEDWRITERDICT) {
		defer file.Close()
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
	defer f.Close()
	_, err = bgzf.HasEOF(f)
	check(err)

	bamReader, err := bam.NewReader(f, THREADNB)
	defer bamReader.Close()
	check(err)

	header := bamReader.Header()
	loadCellIDIndexAndBAMWriter(NUCLEIINDEX, header)

	for _, file := range(WRITERDICT) {
		defer file.Close()
	}

	for _, file := range(BAMWRITERDICT) {
		defer file.Close()
	}

	count := 0

	for {
		record, err = bamReader.Read()
		count++

		switch err {
		case io.EOF:
			break
		case nil:
		default:
			fmt.Printf("ERROR: %s\n",err)
			break
		}

		if record == nil {
			break
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

/*DivideBam divide a bed file */
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

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0, false)
	bedWriter := utils.ReturnWriter(FILENAMEOUT, 6, false)

	defer file.Close()
	defer bedWriter.Close()

	fWrite, err := os.Create(FILENAMEOUT)
	check(err)
	defer fWrite.Close()

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
	defer f.Close()
	_, err = bgzf.HasEOF(f)
	check(err)

	fWrite, err := os.Create(FILENAMEOUT)
	check(err)
	defer fWrite.Close()

	bamReader, err := bam.NewReader(f, THREADNB)
	defer bamReader.Close()
	check(err)

	header := bamReader.Header()
	bamWriter, err := bam.NewWriter(fWrite, header, THREADNB)
	check(err)
	defer bamWriter.Close()
	count := 0

	for {
		record, err = bamReader.Read()
		count++

		switch err {
		case io.EOF:
			break
		case nil:
		default:
			fmt.Printf("ERROR: %s\n",err)
			break
		}

		if record == nil {
			break
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
	defer f.Close()
	scanner := bufio.NewScanner(f)

	CELLIDDICT = make(map[string]bool)

	for scanner.Scan() {
		line := scanner.Text()

		CELLIDDICT[line] = true
	}
}

func loadCellIDIndexAndBAMWriter(fname string, header *sam.Header) {
	f, err := os.Open(fname)
	check(err)
	defer f.Close()
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
	defer f.Close()
	scanner := bufio.NewScanner(f)

	CELLIDDICTMULTIPLE = make(map[string]map[string]bool)
	BEDWRITERDICT = make(map[string]io.WriteCloser)
	OUTFILENAMELIST = make(map[string]bool)

	for scanner.Scan() {
		line := scanner.Text()
		split := strings.Split(line, "\t")
		cellid := split[0]
		filename := split[1]

		ext := filename[len(filename)-7:len(filename)]

		if ext != ".bed.gz" {
			filePath = fmt.Sprintf("%s/%s.bed.gz", OUTPUTDIR, filename)
		} else {
			filePath = fmt.Sprintf("%s/%s", OUTPUTDIR, filename)
		}

		if _, isInside := BEDWRITERDICT[filename]; !isInside {
			BEDWRITERDICT[filename] = utils.ReturnWriter(filePath, 6, false)
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
