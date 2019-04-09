/* Suite of functions dedicated to analyze intersection with genomic regions from a peak file (in bed format) */

package main


import(
	"bufio"
	"log"
	"flag"
	"os"
	"path"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"github.com/biogo/store/interval"
	"fmt"
	"strings"
	"strconv"
	"time"
	"bytes"
)


/*BEDFILENAME bam file name (input) */
var BEDFILENAME string

/*PEAKFILE bam file name (input) */
var PEAKFILE string

/*CELLSIDFNAME file name file with ordered cell IDs (one ID per line) */
var CELLSIDFNAME string

/*CREATESPARSEMATRIX create sparse matrix */
var CREATESPARSEMATRIX bool

/*CELLIDDICT cell ID<->pos */
var CELLIDDICT map[string]uint

/*PEAKIDDICT peak ID<->pos */
var PEAKIDDICT map[string]uint

/*CHRINTERVALDICT peak ID<->pos */
var CHRINTERVALDICT map[string]*interval.IntTree

/*CHRINTERVALDICT peak ID<->pos */
var INTERVALMAPPING map[uintptr]string

/*CELLFEATURE peak ID<->pos */
var CELLFEATURE map[string]*interval.IntTree

/*BOOLSPARSEMATRIX cell x feature sparse matrix  */
var BOOLSPARSEMATRIX map[uint]map[uint]bool

/*FILENAMEOUT  output file name output */
var FILENAMEOUT string


//IntInterval Integer-specific intervals
type IntInterval struct {
	Start, End int
	UID        uintptr
	Payload    interface{}
}


//Overlap rule for two Interval
func (i IntInterval) Overlap(b interval.IntRange) bool {
	// Search for intersection
	return i.End >= b.Start && i.Start <= b.End
}

//ID Return the ID of Interval
func (i IntInterval) ID() uintptr {
	return i.UID
}
//Range Return the range of Interval
func (i IntInterval) Range() interval.IntRange {
	return interval.IntRange{i.Start, i.End}
}

//String Return the string re[ of Interval
func (i IntInterval) String() string {
	return fmt.Sprintf("(%d, %d) id: %d ####\n", i.Start, i.End, i.ID())
}


func main() {
	flag.StringVar(&FILENAMEOUT, "out", "", "name of the output file")
	flag.StringVar(&BEDFILENAME, "bed", "", "name of the bed file")
	flag.StringVar(&PEAKFILE, "peak", "", "name of the bed file")
	flag.StringVar(&CELLSIDFNAME, "cellsID", "", "name of the file containing the ordered list of cell IDs (one ID per line)")
	flag.BoolVar(&CREATESPARSEMATRIX, "coo", false,
		`transform one (-bed) or multiple (use multiple -beds option) into a boolean sparse matrix in COO format
                USAGE: MatUtils -coo -bed  <bedFile> -peaks <bedFile> -cellsID <fname>`)
	flag.Parse()

	if FILENAMEOUT == "" {
		FILENAMEOUT = fmt.Sprintf("%s.tsv", BEDFILENAME)
	}

	tStart := time.Now()

	switch {
	case CREATESPARSEMATRIX:
		switch {
		case BEDFILENAME == "":
			log.Fatal("Error at least one bed file must be provided!")
		case PEAKFILE == "":
			log.Fatal("Error peak file (bed format) must be provided!")
		case CELLSIDFNAME == "":
			log.Fatal("Error cellsID file must be provided!")
		default:
			createBoolSparseMatrix()
		}
	}

	tDiff := time.Now().Sub(tStart)
	fmt.Printf("done in time: %f s \n", tDiff.Seconds())
}


func createBoolSparseMatrix(){
	loadCellIDDict(CELLSIDFNAME)
	loadPeaks(PEAKFILE)

	fmt.Printf("Creating peak interval index...\n")
	tStart := time.Now()
	createPeakIntervalTree()
	tDiff := time.Now().Sub(tStart)
	fmt.Printf("Create peak index done in time: %f s \n", tDiff.Seconds())

	BOOLSPARSEMATRIX = make(map[uint]map[uint]bool)

	for _, pos := range CELLIDDICT {
		BOOLSPARSEMATRIX[pos] = make(map[uint]bool)
	}

	createBoolSparseMatrixOneFile(BEDFILENAME)
	writeBoolMatrixToFile(FILENAMEOUT)
}

func writeBoolMatrixToFile(outfile string) {
	var buffer bytes.Buffer
	var cellPos, featPos uint

	writer, err := os.Create(outfile)
	defer writer.Close()
	check(err)

	for cellPos = range BOOLSPARSEMATRIX {
		for featPos = range BOOLSPARSEMATRIX[cellPos] {
			buffer.WriteString(strconv.Itoa(int(cellPos)))
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.Itoa(int(featPos)))
			buffer.WriteRune('\t')
			buffer.WriteRune('1')
			buffer.WriteRune('\n')

			writer.Write(buffer.Bytes())
			buffer.Reset()
		}
	}

	fmt.Printf("file: %s created!\n", outfile)
}


func createBoolSparseMatrixOneFile(bedfilename string) {
	var line string
	var split []string
	var isInside bool
	var start, end int
	var err error
	var nbReads uint
	var intervals []interval.IntInterface
	var interval interval.IntInterface
	var cellPos,featPos uint

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0, false)

	defer file.Close()

	for bedReader.Scan() {
		line = bedReader.Text()
		nbReads++

		split = strings.Split(line, "\t")

		if cellPos, isInside = CELLIDDICT[split[3]];!isInside {
			continue
		}

		start, err = strconv.Atoi(split[1])
		check(err)

		end, err = strconv.Atoi(split[2])
		check(err)

		if _, isInside = CHRINTERVALDICT[split[0]];!isInside {
			continue
		}

		intervals = CHRINTERVALDICT[split[0]].Get(IntInterval{Start: start, End: end})

		for _, interval = range intervals {
			featPos = PEAKIDDICT[INTERVALMAPPING[interval.ID()]]

			BOOLSPARSEMATRIX[cellPos][featPos] = true
		}
	}
}




func loadCellIDDict(fname string) {
	f, err := os.Open(fname)
	check(err)
	defer f.Close()
	scanner := bufio.NewScanner(f)
	var count uint

	CELLIDDICT = make(map[string]uint)
	count = 0

	for scanner.Scan() {
		line := scanner.Text()
		CELLIDDICT[line] = count
		count++
	}
}

func loadPeaks(fname string) {
	var scanner *bufio.Scanner
	var f *os.File
	var err error

	ext := path.Ext(fname)

	switch ext {
	case ".gz":
		scanner, f = utils.ReturnReader(fname, 0, false)
	default:
		f, err = os.Open(fname)
		check(err)
		scanner = bufio.NewScanner(f)
	}
	defer f.Close()


	var count uint

	PEAKIDDICT = make(map[string]uint)
	count = 0

	for scanner.Scan() {
		line := scanner.Text()
		PEAKIDDICT[line] = count
		count++
	}
}

func createPeakIntervalTree() {
	var split []string
	var chroStr string

	var start, end int
	var err error
	var isInside bool

	CHRINTERVALDICT = make(map[string]*interval.IntTree)
	INTERVALMAPPING = make(map[uintptr]string)

	for key, pos := range PEAKIDDICT {
		split = strings.Split(key, "\t")
		chroStr = split[0]

		start, err = strconv.Atoi(split[1])
		check(err)

		end, err = strconv.Atoi(split[2])
		check(err)

		int := IntInterval{
			Start: start, End: end}
		int.UID = uintptr(uintptr(pos))

		if _, isInside = CHRINTERVALDICT[chroStr];!isInside {
			CHRINTERVALDICT[chroStr] = &interval.IntTree{}
		}

		err = CHRINTERVALDICT[chroStr].Insert(int, false)
		check(err)

		INTERVALMAPPING[int.ID()] = key
	}
}

func check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
