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
)


/*BAMFILENAME bam file name (input) */
var BAMFILENAME string

/*BAMOUTPUT bam file name (output) */
var BAMOUTPUT string

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

/*CELLIDDICTMULTIPLE cell ID<->dict */
var CELLIDDICTMULTIPLE map[string]map[string]bool

/*WRITERDICT filename<->dict */
var WRITERDICT map[string]*os.File

/*BAMWRITERDICT filename<->dict */
var BAMWRITERDICT map[string]*bam.Writer

/*DIVIDE  dividing the bam file tool */
var DIVIDE bool


func main() {
	flag.StringVar(&BAMFILENAME, "bam", "", "name of the bam file")
	flag.StringVar(&BAMOUTPUT, "out", "", "name of the bam file")
	flag.StringVar(&NUCLEIFILE, "cellsID", "", "file with cell IDs")
	flag.StringVar(&OUTPUTDIR, "output_dir", "", "output directory")
	flag.StringVar(&NUCLEIINDEX, "cell_index", "", "nuclei <-> output files index")
	flag.BoolVar(&DIVIDE, "divide", false, "divide the bam file according to barcode file list")
	flag.IntVar(&THREADNB, "thread", 4, "threads concurrency for reading bam file")
	flag.Parse()

	if OUTPUTDIR != "" {
		BAMOUTPUT = fmt.Sprintf("%s/%s", OUTPUTDIR, BAMOUTPUT)
	}

	if BAMFILENAME == "" {
		panic("-bam must be specified!")
	}

	switch {
	case DIVIDE && NUCLEIINDEX!="":
		DivideMultiple()
	case DIVIDE:
		Divide()
	}
}

/*DivideMultiple divide the bam file */
func DivideMultiple() {
	var record * sam.Record
	var read string
	var readID string
	var bamWriter *bam.Writer
	var isInside bool
	var filename string

	tStart := time.Now()

	f, err := os.Open(BAMFILENAME)
	check(err)
	defer f.Close()
	_, err = bgzf.HasEOF(f)
	check(err)

	bamReader, err := bam.NewReader(f, THREADNB)
	defer bamReader.Close()
	check(err)

	header := bamReader.Header()
	loadCellIDIndex(NUCLEIINDEX, header)
	check(err)

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

	tDiff := time.Now().Sub(tStart)
	fmt.Printf("done in time: %f s for %d reads\n", tDiff.Seconds(), count-1)
}


/*Divide divide the bam file */
func Divide() {
	var record * sam.Record
	var read string
	var readID string

	if BAMOUTPUT == "" {
		panic("-out must be specified!")
	}

	if NUCLEIFILE == "" {
		panic("-cellsID must be specified!")
	}

	tStart := time.Now()

	loadCellIDDict(NUCLEIFILE)

	f, err := os.Open(BAMFILENAME)
	check(err)
	defer f.Close()
	_, err = bgzf.HasEOF(f)
	check(err)

	fWrite, err := os.Create(BAMOUTPUT)
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

	tDiff := time.Now().Sub(tStart)
	fmt.Printf("done in time: %f s for %d reads\n", tDiff.Seconds(), count-1)
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

func loadCellIDIndex(fname string, header *sam.Header) {
	f, err := os.Open(fname)
	check(err)
	defer f.Close()
	scanner := bufio.NewScanner(f)
	var filePath string

	CELLIDDICTMULTIPLE = make(map[string]map[string]bool)
	WRITERDICT = make(map[string]*os.File)
	BAMWRITERDICT = make(map[string]*bam.Writer)

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
