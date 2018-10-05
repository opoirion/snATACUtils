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
)


/*BAMFILENAME bam file name (input) */
var BAMFILENAME string

/*BAMOUTPUT bam file name (output) */
var BAMOUTPUT string

/*THREADNB number of threads for reading the bam file */
var THREADNB int

/*NUCLEIFILE file with cells ID */
var NUCLEIFILE string

/*CELLIDDICT cell ID<->dict */
var CELLIDDICT map[string]bool

/*DIVIDE  dividing the bam file tool */
var DIVIDE bool


func main() {
	flag.StringVar(&BAMFILENAME, "bam", "", "name of the bam file")
	flag.StringVar(&BAMOUTPUT, "out", "", "name of the bam file")
	flag.StringVar(&NUCLEIFILE, "cellsID", "", "file with cell IDs")
	flag.BoolVar(&DIVIDE, "divide", false, "dividing the bam files")
	flag.IntVar(&THREADNB, "thread", 4, "threads concurrency for reading bam file")
	flag.Parse()


	if BAMFILENAME == "" {
		panic("-bam must be specified!")
	}

	switch {
	case DIVIDE:
		Divide()
	}
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

func check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
