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
	"sync"
)


/*INFILES multiple input files */
var INFILES utils.ArrayFlags

/*BEDFILENAME bed file name (input) */
var BEDFILENAME string

/*PEAKFILE bed file name (input) */
var PEAKFILE string

/*CELLSIDFNAME file name file with ordered cell IDs (one ID per line) */
var CELLSIDFNAME string

/*CREATECOOMATRIX create COO sparse matrix */
var CREATECOOMATRIX bool

/*READINPEAK read in peak */
var READINPEAK bool

/*MERGEOUTPUTS merge output files */
var MERGEOUTPUTS bool

/*SEP separator for writing output */
var SEP string

/*CELLIDCOUNT cell ID<->count */
var CELLIDCOUNT map[string]int

/*CELLIDDICT cell ID<->pos */
var CELLIDDICT map[string]uint

/*MUTEX global mutex */
var MUTEX *sync.Mutex

/*CELLMUTEXDICT feature pos<->sync.Mutex */
var CELLMUTEXDICT map[uint]*sync.Mutex

/*BOOLSPARSEMATRIX cell x feature sparse matrix  */
var BOOLSPARSEMATRIX map[uint]map[uint]bool

/*FILENAMEOUT  output file name output */
var FILENAMEOUT string

/*THREADNB number of threads for reading the bam file */
var THREADNB int

/*BUFFERSIZE buffer size for multithreading */
const BUFFERSIZE = 1000000


func main() {
	flag.StringVar(&FILENAMEOUT, "out", "", "name of the output file")
	flag.StringVar(&BEDFILENAME, "bed", "", "name of the bed file")
	flag.Var(&INFILES, "in", "name of the input file(s)")
	flag.StringVar(&PEAKFILE, "ygi", "", "name of the bed file containing the region of interest( i.e. PEAK )")
	flag.StringVar(&CELLSIDFNAME, "xgi", "", "name of the file containing the ordered list of cell IDs (one ID per line)")
	flag.StringVar(&SEP, "delimiter", "\t", "delimiter used to write the output file (default \t)")
	flag.BoolVar(&CREATECOOMATRIX, "coo", false,
		`transform one (-bed) or multiple (use multiple -beds option) into a boolean sparse matrix in COO format
                USAGE: ATACMatTools -coo -bed  <bedFile> -ygi <bedFile> -xgi <fname>`)
	flag.BoolVar(&MERGEOUTPUTS, "merge", false,
		`merge multiple matrices results into one output file
                USAGE: ATACMatTools -coo -merge -xgi <fname> -in <matrixFile1> -in <matrixFile2> ...`)
	flag.BoolVar(&READINPEAK, "count", false,
		`Count the number of reads in peaks for each cell
                USAGE: ATACMatTools -count  -xgi <fname> -ygi <bedfile> -bed <bedFile>`)
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
	case CREATECOOMATRIX || READINPEAK:
		switch {
		case CELLSIDFNAME == "" && !READINPEAK:
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
		case READINPEAK:
			computeReadsInPeaksForCell()
		default:
			createBoolSparseMatrix()
		}
	}

	tDiff := time.Since(tStart)
	fmt.Printf("done in time: %f s \n", tDiff.Seconds())
}


func computeReadsInPeaksForCell(){
	fmt.Printf("load indexes...\n")

	if CELLSIDFNAME != "" {
		loadCellIDDict(CELLSIDFNAME)
	}

	utils.LoadPeaks(PEAKFILE)
	utils.CreatePeakIntervalTree()

	CELLIDCOUNT = make(map[string]int)

	fmt.Printf("init mutexes...\n")
	initMutexDict()
	fmt.Printf("init threading...\n")
	utils.InitIntervalDictsThreading(THREADNB)
	createReadInPeakOneFileThreading(BEDFILENAME)
	writeCellCounter(FILENAMEOUT)
}


func initBoolSparseMatrix() {
	BOOLSPARSEMATRIX = make(map[uint]map[uint]bool)

	for _, pos := range CELLIDDICT {
		BOOLSPARSEMATRIX[pos] = make(map[uint]bool)
	}
}


func createBoolSparseMatrix(){
	fmt.Printf("load indexes...\n")
	loadCellIDDict(CELLSIDFNAME)
	utils.LoadPeaks(PEAKFILE)

	initBoolSparseMatrix()
	utils.CreatePeakIntervalTree()

	switch{
	case THREADNB > 1:
		fmt.Printf("init mutexes...\n")
		initMutexDict()
		fmt.Printf("init threading...\n")
		utils.InitIntervalDictsThreading(THREADNB)
		fmt.Printf("launching sparse matrices creation...\n")
		createBoolSparseMatrixOneFileThreading(BEDFILENAME)
	default:
		createBoolSparseMatrixOneFile(BEDFILENAME)
	}

	writeBoolMatrixToFile(FILENAMEOUT)
}

func initMutexDict() {
	MUTEX = &sync.Mutex{}

	CELLMUTEXDICT = make(map[uint]*sync.Mutex)

	for _, pos := range CELLIDDICT {
		CELLMUTEXDICT[pos] = &sync.Mutex{}
	}
}


func writeBoolMatrixToFile(outfile string) {
	fmt.Printf("writing to output file...\n")

	var buffer bytes.Buffer
	var cellPos, featPos uint

	writer := utils.ReturnWriter(outfile)

	defer utils.CloseFile(writer)

	for cellPos = range BOOLSPARSEMATRIX {
		for featPos = range BOOLSPARSEMATRIX[cellPos] {
			buffer.WriteString(strconv.Itoa(int(cellPos)))
			buffer.WriteString(SEP)
			buffer.WriteString(strconv.Itoa(int(featPos)))
			buffer.WriteString(SEP)
			buffer.WriteRune('1')
			buffer.WriteRune('\n')

			writer.Write(buffer.Bytes())
			buffer.Reset()
		}
	}

	fmt.Printf("file: %s created!\n", outfile)
}

func writeCellCounter(outfile string) {
	var buffer bytes.Buffer
	var cellID string

	writer := utils.ReturnWriter(outfile)

	defer utils.CloseFile(writer)

	for cellID = range CELLIDCOUNT {
		buffer.WriteString(cellID)
		buffer.WriteString(SEP)
		buffer.WriteString(strconv.Itoa(int(CELLIDCOUNT[cellID])))
		buffer.WriteRune('\n')

		writer.Write(buffer.Bytes())
		buffer.Reset()
	}

	fmt.Printf("file: %s created!\n", outfile)
}

/*mergeCOOFiles merge multiple COO output files*/
func mergeCOOFiles(filenames []string) {
	fmt.Printf("creating xgi index..\n")
	loadCellIDDict(CELLSIDFNAME)
	initBoolSparseMatrix()

	for _, filename := range filenames {
		fmt.Printf("merging file: %s\n", filename)
		mergeCOOFile(filename)
	}

	writeBoolMatrixToFile(FILENAMEOUT)
}


/*mergeCOOFile add one file to the matrix*/
func mergeCOOFile(filename string) {
	var scanner *bufio.Scanner
	var f *os.File
	var err error
	var split []string
	var start, stop int

	scanner, f = utils.ReturnReader(filename, 0)

	defer utils.CloseFile(f)

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), SEP)
		start, err = strconv.Atoi(split[0])
		utils.Check(err)
		stop, err = strconv.Atoi(split[1])
		utils.Check(err)

		BOOLSPARSEMATRIX[uint(start)][uint(stop)] = true
	}
}


/*createBoolSparseMatrixOneFileThreading ceate the bool Sparse Matrix for one bed file using multi-threading*/
func createBoolSparseMatrixOneFileThreading(bedfilename string) {
	var nbReads uint
	var bufferLine1 [BUFFERSIZE]string
	var bufferLine2 [BUFFERSIZE]string
	var bufferPointer * [BUFFERSIZE]string

	isBuffer1 := true
	bufferPointer = &bufferLine1

	var bufferIt int
	var waiting sync.WaitGroup

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0)

	defer utils.CloseFile(file)

	scanBed:
	for bedReader.Scan() {
		bufferPointer[bufferIt] = bedReader.Text()
		nbReads++
		bufferIt++


		if bufferIt >= BUFFERSIZE {
			chunk := bufferIt / THREADNB
			bufferStart := 0
			bufferStop := chunk

			for i := 0; i < THREADNB;i++{

				waiting.Add(1)
				go updateBoolSparseMatrixOneThread(bufferPointer , bufferStart, bufferStop, i, &waiting)

				bufferStart += chunk
				bufferStop += chunk

				if i == THREADNB - 1 {
					bufferStop = bufferIt
				}
			}

			bufferIt = 0

			if isBuffer1 {
				bufferPointer = &bufferLine2
				isBuffer1 = false
				goto scanBed
			} else {
				bufferPointer = &bufferLine1
				isBuffer1 = true
				waiting.Wait()
			}
		}
	}

	waiting.Wait()

	if bufferIt > 0 {
		waiting.Add(1)
		updateBoolSparseMatrixOneThread(bufferPointer , 0, bufferIt, 0, &waiting)
	}
}

func createReadInPeakOneFileThreading(bedfilename string) {
	var nbReads uint
	var bufferLine1 [BUFFERSIZE]string
	var bufferLine2 [BUFFERSIZE]string
	var bufferPointer * [BUFFERSIZE]string

	var bufferIt int
	var waiting sync.WaitGroup

	isBuffer1 := true
	bufferPointer = &bufferLine1

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0)

	defer utils.CloseFile(file)

	for bedReader.Scan() {
		bufferPointer[bufferIt] = bedReader.Text()
		nbReads++
		bufferIt++

		scanBed:
		if bufferIt >= BUFFERSIZE {
			chunk := bufferIt / THREADNB
			bufferStart := 0
			bufferStop := chunk

			for i := 0; i < THREADNB;i++{

				waiting.Add(1)
				go updateReadInPeakThread(bufferPointer , bufferStart, bufferStop, i, &waiting)

				bufferStart += chunk
				bufferStop += chunk

				if i == THREADNB - 1 {
					bufferStop = bufferIt
				}
			}

			bufferIt = 0

			if isBuffer1 {
				bufferPointer = &bufferLine2
				isBuffer1 = false
				goto scanBed
			} else {
				bufferPointer = &bufferLine1
				isBuffer1 = true
				waiting.Wait()
			}
		}
	}

	waiting.Wait()

	if bufferIt > 0 {
		waiting.Add(1)
		updateReadInPeakThread(bufferPointer , 0, bufferIt, 0, &waiting)
	}
}

func updateReadInPeakThread(bufferLine * [BUFFERSIZE]string, bufferStart ,bufferStop, threadnb int,
	waiting * sync.WaitGroup) {
	defer waiting.Done()
	var split []string
	var isInside bool
	var start, end int
	var err error

	isCellsID := true

	if CELLSIDFNAME == "" {
		isCellsID = false
	}

	var intervals []interval.IntInterface

	for i := bufferStart; i < bufferStop;i++ {

		split = strings.Split(bufferLine[i], "\t")

		if isCellsID {

			if _, isInside = CELLIDDICT[split[3]];!isInside {
				continue
			}

		}

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		if _, isInside = utils.CHRINTERVALDICT[split[0]];!isInside {
			continue
		}

		intervals = utils.CHRINTERVALDICTTHREAD[threadnb][split[0]].Get(
			utils.IntInterval{Start: start, End: end})

		MUTEX.Lock()

		for range intervals {
			CELLIDCOUNT[split[3]]++
		}

		MUTEX.Unlock()
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
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		if _, isInside = utils.CHRINTERVALDICT[split[0]];!isInside {
			continue
		}

		intervals = utils.CHRINTERVALDICTTHREAD[threadnb][split[0]].Get(
			utils.IntInterval{Start: start, End: end})

		CELLMUTEXDICT[cellPos].Lock()

		for _, int = range intervals {
			featPos = utils.PEAKIDDICT[utils.INTERVALMAPPING[int.ID()]]
			BOOLSPARSEMATRIX[cellPos][featPos] = true
		}

		CELLMUTEXDICT[cellPos].Unlock()
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

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0)

	defer utils.CloseFile(file)

	for bedReader.Scan() {
		line = bedReader.Text()
		nbReads++

		split = strings.Split(line, "\t")

		if cellPos, isInside = CELLIDDICT[split[3]];!isInside {
			continue
		}

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		if _, isInside = utils.CHRINTERVALDICT[split[0]];!isInside {
			continue
		}

		intervals = utils.CHRINTERVALDICT[split[0]].Get(
			utils.IntInterval{Start: start, End: end})

		for _, interval = range intervals {
			featPos = utils.PEAKIDDICT[utils.INTERVALMAPPING[interval.ID()]]

			BOOLSPARSEMATRIX[cellPos][featPos] = true
		}
	}
}



/*loadCellIDDict ...*/
func loadCellIDDict(fname string) {
	f, err := os.Open(fname)
	utils.Check(err)
	defer utils.CloseFile(f)
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
