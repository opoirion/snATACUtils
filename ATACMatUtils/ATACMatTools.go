/* Suite of functions dedicated to analyze intersection with genomic regions from a peak file (in bed format) */

package main


import(
	"bufio"
	"log"
	"flag"
	"os"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"github.com/biogo/store/interval"
	"fmt"
	"strings"
	"strconv"
	"time"
	"bytes"
	"io"
	"sync"
	"github.com/jinzhu/copier"
)

/*INFILES multiple input files */
var INFILES utils.ArrayFlags

/*BEDFILENAME bam file name (input) */
var BEDFILENAME string

/*PEAKFILE bam file name (input) */
var PEAKFILE string

/*CELLSIDFNAME file name file with ordered cell IDs (one ID per line) */
var CELLSIDFNAME string

/*CREATECOOMATRIX create COO sparse matrix */
var CREATECOOMATRIX bool

/*MERGECOOOUTPUTS merge output files */
var MERGEOUTPUTS bool

/*CELLIDDICT cell ID<->pos */
var CELLIDDICT map[string]uint

/*PEAKIDDICT peak ID<->pos */
var PEAKIDDICT map[string]uint

/*CHRINTERVALDICT peak ID<->pos */
var CHRINTERVALDICT map[string]*interval.IntTree

/*CHRINTERVALDICTTHREAD peak ID<->pos */
var CHRINTERVALDICTTHREAD map[int]map[string]*interval.IntTree

/*FEATMUTEXDICT feature pos<->sync.Mutex */
var FEATMUTEXDICT map[uint]*sync.Mutex

/*CHRINTERVALDICT peak ID<->pos */
var INTERVALMAPPING map[uintptr]string

/*BOOLSPARSEMATRIX cell x feature sparse matrix  */
var BOOLSPARSEMATRIX map[uint]map[uint]bool

/*FILENAMEOUT  output file name output */
var FILENAMEOUT string

/*THREADNB number of threads for reading the bam file */
var THREADNB int

/*BUFFERSIZE buffer size for multithreading */
const BUFFERSIZE = 1000000

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
	flag.Var(&INFILES, "in", "name of the input file(s)")
	flag.StringVar(&PEAKFILE, "ygi", "", "name of the bed file containing the region of interest( i.e. PEAK )")
	flag.StringVar(&CELLSIDFNAME, "xgi", "", "name of the file containing the ordered list of cell IDs (one ID per line)")
	flag.BoolVar(&CREATECOOMATRIX, "coo", false,
		`transform one (-bed) or multiple (use multiple -beds option) into a boolean sparse matrix in COO format
                USAGE: ATACMatTools -coo -bed  <bedFile> -ygi <bedFile> -xgi <fname>`)
	flag.BoolVar(&MERGEOUTPUTS, "merge", false,
		`merge multiple matrices results into one output file
                USAGE: ATACMatTools -coo -merge -xgi <fname> -in <matrixFile1> -in <matrixFile2> ...`)
	flag.IntVar(&THREADNB, "threads", 1, "threads concurrency")
	flag.Parse()

	switch {
	case FILENAMEOUT == "" && len(INFILES) > 0:
		FILENAMEOUT = fmt.Sprintf("%s.coo.gz", INFILES[0])
	case FILENAMEOUT == "" && BEDFILENAME != "":
		FILENAMEOUT = fmt.Sprintf("%s.coo.gz", BEDFILENAME)
	case FILENAMEOUT == "":
		FILENAMEOUT = fmt.Sprintf("output.coo.gz")
	}

	tStart := time.Now()

	switch {
	case CREATECOOMATRIX:
		switch {
		case CELLSIDFNAME == "":
			log.Fatal("Error -xgi file must be provided!")
		case MERGEOUTPUTS:
			if len(INFILES) == 0 {
				log.Fatal("Error at least one input (-in) file must be provided!")
			}


			mergeCOOFiles(INFILES)
		case BEDFILENAME == "":
			log.Fatal("Error at least one bed file must be provided!")
		case PEAKFILE == "":
			log.Fatal("Error peak file -ygi (bed format) must be provided!")
		default:
			createBoolSparseMatrix()
		}
	}

	tDiff := time.Now().Sub(tStart)
	fmt.Printf("done in time: %f s \n", tDiff.Seconds())
}


func createBoolSparseMatrix(){
	fmt.Printf("load indexes...\n")
	loadCellIDDict(CELLSIDFNAME)
	loadPeaks(PEAKFILE)

	BOOLSPARSEMATRIX = make(map[uint]map[uint]bool)

	for _, pos := range CELLIDDICT {
		BOOLSPARSEMATRIX[pos] = make(map[uint]bool)
	}
	fmt.Printf("create peak interval tree...\n")
	tStart := time.Now()
	createPeakIntervalTree()
	tDiff := time.Now().Sub(tStart)
	fmt.Printf("Create peak index done in time: %f s \n", tDiff.Seconds())

	switch{
	case THREADNB > 1:
		fmt.Printf("init mutexes...\n")
		initMutexDict()
		fmt.Printf("init threading...\n")
		initIntervalDictsThreading()
		fmt.Printf("launching sparse matrices creation...\n")
		createBoolSparseMatrixOneFileThreading(BEDFILENAME)
	default:
		createBoolSparseMatrixOneFile(BEDFILENAME)
	}

	fmt.Printf("writing to output file...\n")
	writeBoolMatrixToFile(FILENAMEOUT)
}

func initMutexDict() {
	FEATMUTEXDICT = make(map[uint]*sync.Mutex)

	for _, pos := range PEAKIDDICT {
		FEATMUTEXDICT[pos] = &sync.Mutex{}
	}
}

func initIntervalDictsThreading() {
	CHRINTERVALDICTTHREAD = make(map[int]map[string]*interval.IntTree)

	for i := 0;i< THREADNB;i++ {
		CHRINTERVALDICTTHREAD[i] = make(map[string]*interval.IntTree)

		for key, tree := range CHRINTERVALDICT {
			CHRINTERVALDICTTHREAD[i][key] = &interval.IntTree{}
			copier.Copy(CHRINTERVALDICTTHREAD[i][key], tree)
		}
	}
}

func writeBoolMatrixToFile(outfile string) {
	var buffer bytes.Buffer
	var cellPos, featPos uint
	var writer io.WriteCloser

	writer = utils.ReturnWriter(outfile, 0, false)

	defer writer.Close()

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

/*mergeCOOFiles merge multiple COO output files*/
func mergeCOOFiles(filenames []string) {
	fmt.Printf("creating xgi index..\n")
	loadCellIDDict(CELLSIDFNAME)

	BOOLSPARSEMATRIX = make(map[uint]map[uint]bool)

	for _, pos := range CELLIDDICT {
		BOOLSPARSEMATRIX[pos] = make(map[uint]bool)
	}

	for _, filename := range filenames {
		fmt.Printf("merging file: %s\n", filename)
		mergeCOOFile(filename)
	}

	fmt.Printf("writing to output file...\n")

	writeBoolMatrixToFile(FILENAMEOUT)
}


/*mergeCOOFile add one file to the matrix*/
func mergeCOOFile(filename string) {
	var scanner *bufio.Scanner
	var f *os.File
	var err error
	var split []string
	var start, stop int

	scanner, f = utils.ReturnReader(filename, 0, false)

	defer f.Close()

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), "\t")
		start, err = strconv.Atoi(split[0])
		check(err)
		stop, err = strconv.Atoi(split[1])
		check(err)

		BOOLSPARSEMATRIX[uint(start)][uint(stop)] = true
	}
}


/*createBoolSparseMatrixOneFileThreading ceate the bool Sparse Matrix for one bed file using multi-threading*/
func createBoolSparseMatrixOneFileThreading(bedfilename string) {
	var nbReads uint
	var bufferLine [BUFFERSIZE]string
	var bufferIt int
	var waiting sync.WaitGroup

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0, false)

	defer file.Close()

	for bedReader.Scan() {
		bufferLine[bufferIt] = bedReader.Text()
		nbReads++
		bufferIt++

		if bufferIt >= BUFFERSIZE {
			chunk := bufferIt / THREADNB
			bufferStart := 0
			bufferStop := chunk

			for i := 0; i < THREADNB;i++{

				waiting.Add(1)
				go updateBoolSparseMatrixOneThread(&bufferLine , bufferStart, bufferStop, i, &waiting)

				bufferStart += chunk
				bufferStop += chunk

				if i == THREADNB - 1 {
					bufferStop = bufferIt
				}
			}

			bufferIt = 0
		}

		waiting.Wait()
	}

	if bufferIt > 0 {
		waiting.Add(1)
		updateBoolSparseMatrixOneThread(&bufferLine , 0, bufferIt, 0, &waiting)
	}
}


func updateBoolSparseMatrixOneThread(bufferLine * [BUFFERSIZE]string, bufferStart ,bufferStop, threadnb int,
	waiting * sync.WaitGroup) {
	defer waiting.Done()
	var split []string
	var isInside bool
	var start, end int
	var err error

	var intervals []interval.IntInterface
	var int interval.IntInterface
	var cellPos,featPos uint

	for i := bufferStart; i < bufferStop;i++ {

		split = strings.Split(bufferLine[i], "\t")

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

		intervals = CHRINTERVALDICTTHREAD[threadnb][split[0]].Get(IntInterval{Start: start, End: end})

		for _, int = range intervals {
			featPos = PEAKIDDICT[INTERVALMAPPING[int.ID()]]
			FEATMUTEXDICT[featPos].Lock()
			BOOLSPARSEMATRIX[cellPos][featPos] = true
			FEATMUTEXDICT[featPos].Unlock()
		}
	}
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



/*loadCellIDDict ...*/
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

/*loadPeaks ...*/
func loadPeaks(fname string) {
	var scanner *bufio.Scanner
	var f *os.File

	scanner, f = utils.ReturnReader(fname, 0, false)

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

/*createPeakIntervalTree ...*/
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
