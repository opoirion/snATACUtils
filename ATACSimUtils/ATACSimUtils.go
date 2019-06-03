/* Suite of functions dedicated to generate Simulated snATAC-Seq data */

package main


import(
	"log"
	"flag"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"fmt"
	"github.com/valyala/fastrand"
	"sync"
	"path"
	"io"
	"bytes"
	"strings"
	"strconv"
	"time"
	"math/rand"
)


/*BEDFILENAME bed file name (input) */
var BEDFILENAME string

/*BEDFILENAMES multiple input files */
var BEDFILENAMES utils.ArrayFlags

/*FILENAMEOUT  output file name output */
var FILENAMEOUT string

/*TAGNAME  tag name for the simulated cells */
var TAGNAME string

/*MEAN  mean of the nb of reads per cell dist */
var MEAN float64

/*STD  std of the nb of reads per cell dist */
var STD float64

/*CELLNB Number of cells to generate */
var CELLNB int

/*THREADNB number of threads for reading the bam file */
var THREADNB int

/*SIMULATEBED  Simulate a scATAC-Seq bed file using a reference bed file*/
var SIMULATEBED bool

/*SEED  Seed used for random processes*/
var SEED int

/*READSARRAY  cell->nb reads array*/
var READSARRAY []float64

/*BUFFERARRAY  cell->nb reads array*/
var BUFFERARRAY []bytes.Buffer

/*THREADSCHANNEL  cell->nb reads array*/
var THREADSCHANNEL chan int

/*MUTEX  global mutex*/
var MUTEX sync.Mutex

/*WAITING  waiting group*/
var WAITING * sync.WaitGroup

/*BUFFERSIZE  buffer size*/
const BUFFERSIZE = 50000

/*BUFFERLINEARRAY  cell->nb reads array*/
var BUFFERLINEARRAY [][BUFFERSIZE]string


func main() {
	flag.StringVar(&BEDFILENAME, "bed", "", "name of the bed file")
	flag.StringVar(&TAGNAME, "tag", "", "tag name for the simulated cells")
	flag.Var(&BEDFILENAMES, "beds", "name of several bed files")
	flag.IntVar(&THREADNB, "threads", 1, "threads concurrency")
	flag.IntVar(&SEED, "seed", 2019, "Seed used for random processes")
	flag.IntVar(&CELLNB, "nb", 5000, "Number of cells to generate")
	flag.Float64Var(&MEAN, "mean", 4000, "Average nb. of reads per cell used")
	flag.Float64Var(&STD, "std", 2000, "Std. of the nb. of reads per cell used")
	flag.StringVar(&FILENAMEOUT, "out", "", "name/tag the output file(s)")
	flag.BoolVar(&SIMULATEBED, "simulate", false, `Simulate scATAC-Seq bed files
                      USAGE: ATACSimUtils -simulate -nb <int> -mean <float> std <float> -bed <bedfile> (-threads <int> -out <string> -tag <string>)`)
	flag.Parse()

	BEDFILENAMES = append(BEDFILENAMES, BEDFILENAME)

	switch{
	case SIMULATEBED:
		switch {
		case len(BEDFILENAMES) == 0:
			log.Fatal("Error At least one bed file should be given as input")
		case FILENAMEOUT == "" && len(BEDFILENAMES) > 1:
			FILENAMEOUT = "simulated"
		}

		simulateBedFiles(BEDFILENAMES)
	}
}


func simulateBedFiles(bedfilenames []string) {
	outputfile := FILENAMEOUT
	WAITING = &sync.WaitGroup{}

	for pos, bedfilename := range bedfilenames {
		fmt.Printf("Simulating File nb %d: %s...\n", pos + 1, bedfilename)
		nbLines := countNbLines(bedfilename)
		fmt.Printf("Nb reads: %d\n", nbLines)

		if len(bedfilenames) > 1 {
			ext := path.Ext(bedfilename)
			outputfile = fmt.Sprintf("%s.%s%s",
				FILENAMEOUT[:len(FILENAMEOUT) - len(ext)],
				outputfile, ext)
		}

		simulateOneBedFile(bedfilename, outputfile, nbLines, pos)
		fmt.Printf("File: %s written!\n", outputfile)
	}
}


func simulateOneBedFile(bedfilename, outputfile string, nbLines , it int) {
	tStart := time.Now()

	fmt.Printf("Initating random number of reads per cell...\n")
	initNbReads(it, nbLines)

	BUFFERARRAY = make([]bytes.Buffer, THREADNB)
	BUFFERLINEARRAY = make([][BUFFERSIZE]string, THREADNB)
	THREADSCHANNEL = make(chan int, THREADNB)

	writer := utils.ReturnWriter(outputfile)
	defer writer.Close()

	for i:=0;i<THREADNB;i++ {
		THREADSCHANNEL <- i
	}


	bedReader, file := utils.ReturnReader(bedfilename, 0)
	defer file.Close()
	threadID := <-THREADSCHANNEL
	lineID := 0

	for bedReader.Scan() {
		BUFFERLINEARRAY[threadID][lineID] = bedReader.Text()

		lineID++

		if lineID >= BUFFERSIZE {
			WAITING.Add(1)
			go processOneRead(&BUFFERLINEARRAY[threadID], &writer, threadID, lineID)
			lineID = 0
			threadID = <-THREADSCHANNEL
		}
	}

	WAITING.Add(1)
	go processOneRead(&BUFFERLINEARRAY[threadID], &writer, threadID, lineID)
	WAITING.Wait()

	tDiff := time.Now().Sub(tStart)
	fmt.Printf("Simulating one bed done in time: %f s \n", tDiff.Seconds())
}

func processOneRead(lines * [BUFFERSIZE]string, writer * io.WriteCloser, threadID, lineEnd int) {
	defer WAITING.Done()
	var randNum float64
	var split []string
	var i int

	write := false

	for pos := 0;pos < lineEnd;pos++ {

		split = strings.Split(lines[pos], "\t")

		for i = 0; i < len(READSARRAY); i++ {
			randNum = float64(fastrand.Uint32n(100000)) / 100000.0

			if READSARRAY[i] > randNum {
				write = true

				BUFFERARRAY[threadID].WriteString(split[0])
				BUFFERARRAY[threadID].WriteRune('\t')
				BUFFERARRAY[threadID].WriteString(split[1])
				BUFFERARRAY[threadID].WriteRune('\t')
				BUFFERARRAY[threadID].WriteString(split[2])
				BUFFERARRAY[threadID].WriteRune('\t')
				BUFFERARRAY[threadID].WriteString("SIM")
				BUFFERARRAY[threadID].WriteString(TAGNAME)
				BUFFERARRAY[threadID].WriteString(strconv.Itoa(i))
				BUFFERARRAY[threadID].WriteRune('\n')
			}
		}

	}

	if write {
		MUTEX.Lock()
		(*writer).Write(BUFFERARRAY[threadID].Bytes())
		BUFFERARRAY[threadID].Reset()
		MUTEX.Unlock()
	}

	THREADSCHANNEL <- threadID
}


func initNbReads(it int, nbLines int) {
	rand.Seed(int64(SEED + it))

	READSARRAY = make([]float64, CELLNB, CELLNB)
	nbLinesFloat := float64(nbLines)

	for i :=0;i < len(READSARRAY);i++ {
		READSARRAY[i] = (rand.NormFloat64() * STD + MEAN) / nbLinesFloat
	}
}


func countNbLines(bedfile string) int {
	bedReader, file := utils.ReturnReader(bedfile, 0)

	defer file.Close()
	nbLines := 0

	for bedReader.Scan() {
		nbLines++
	}

	return nbLines
}
