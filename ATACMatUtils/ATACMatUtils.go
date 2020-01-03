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
var BEDFILENAME utils.Filename

/*PEAKFILE bed file name (input) */
var PEAKFILE utils.Filename

/*CELLSIDFNAME file name file with ordered cell IDs (one ID per line) */
var CELLSIDFNAME utils.Filename

/*COO create COO bool sparse matrix */
var COO bool

/*USECOUNT create COO count sparse matrix */
var USECOUNT bool

/*CREATEBINMATRIX create BIN sparse matrix */
var CREATEBINMATRIX bool

/*READINPEAK read in peak */
var READINPEAK bool

/*NORM normalise bin matrix */
var NORM bool

/*MERGEOUTPUTS merge output files */
var MERGEOUTPUTS bool

/*TAIJI Use TAIJI format */
var TAIJI bool

/*SPLIT split computation */
var SPLIT int

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

/*INTSPARSEMATRIX cell x feature sparse matrix  */
var INTSPARSEMATRIX []map[uint]int

/*FILENAMEOUT  output file name output */
var FILENAMEOUT string

/*THREADNB number of threads for reading the bam file */
var THREADNB int

/*YGIOUT number of threads for reading the bam file */
var YGIOUT string

/*YGIDIM dim of ygi when reading a COO file */
var YGIDIM int

/*XGIDIM dim of xgi when reading a COO file */
var XGIDIM int

/*BUFFERSIZE buffer size for multithreading */
const BUFFERSIZE = 1000000


func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `
#################### MODULE TO CREATE (cell x genomic region) SPARSE MATRIX ########################
"""Boolean / interger Peak matrix """
transform one (-bed) or multiple bed files into a sparse matrix
USAGE: ATACMatTools -coo -bed  <bedFile> -ygi <bedFile> -xgi <fname> (-threads <int> -out <fname> -use_count -taiji -bed <fileName2>)

"""Create a cell x bin matrix: -bin """
transform one (-bed) or multiple (use multiple -bed options) bed file into a bin (using float) sparse matrix. If ygi provided, reads intersecting these bin are ignored

USAGE: ATACMatTools -bin -bed  <bedFile> (optional -ygi <bedFile> -xgi <fname> -bin_size <int> -ygi_out <string> -norm -taiji -coo)

"""Count the number of reads in peaks for each cell: -count """
USAGE: ATACMatTools -count  -xgi <fname> -ygi <bedfile> -bed <bedFile>

"""Merge multiple matrices results into one output file: -merge """
USAGE: ATACMatTools -coo/taiji -merge -xgi <fname> -in <matrixFile1> -in <matrixFile2> ... (optional /bin)

`)
		 flag.PrintDefaults()
	}


	flag.IntVar(&SPLIT, "split", 0, "Split computation into n iterative chuncks (to reduce RAM usage for very large matrices)")
	flag.IntVar(&BINSIZE, "bin_size", 5000, "Size of the bin for bin matrix")
	flag.StringVar(&FILENAMEOUT, "out", "", "name of the output file")
	flag.Var(&BEDFILENAME, "bed", "name of the bed file")
	flag.Var(&INFILES, "in", "name of the input file(s)")
	flag.Var(&PEAKFILE, "ygi", "name of the bed file containing the region of interest( i.e. PEAK )")
	flag.Var(&CELLSIDFNAME, "xgi", "name of the file containing the ordered list of cell IDs (one ID per line)")
	flag.StringVar(&SEP, "delimiter", "\t", "delimiter used to write the output file (default \t)")
	flag.StringVar(&YGIOUT, "ygi_out", "", "Write the output bin bed file in this specific file")
	flag.BoolVar(&TAIJI, "taiji", false,
		`Convert COO matrix file to sparse matrix format required for taiji-utils. See https://github.com/Taiji-pipeline/Taiji-utils/`)
	flag.BoolVar(&USECOUNT, "use_count", false,
		`Use read count instead of boolean value`)
	flag.BoolVar(&COO, "coo", false,
		`Use COO format as output`)

	flag.BoolVar(&NORM, "norm", false, "Normalize bin matrix per read depth for each cell")

	flag.BoolVar(&CREATEBINMATRIX, "bin", false,
		`transform one (-bed) or multiple (use multiple -beds option) into a bin (using float) sparse matrix in COO format.`)
	flag.BoolVar(&MERGEOUTPUTS, "merge", false, `merge multiple matrices results into one output file`)
	flag.BoolVar(&READINPEAK, "count", false, `Count the number of reads in peaks for each cell`)
	flag.IntVar(&THREADNB, "threads", 1, "threads concurrency")
	flag.Parse()

	var tag string


	if CREATEBINMATRIX {
		tag = "bin."
	}

	switch {
	case COO:
		tag = fmt.Sprintf("%scoo", tag)
	case TAIJI:
		tag = fmt.Sprintf("%staiji-mat", tag)
	}

	switch {
	case FILENAMEOUT == "" && len(INFILES) > 0:
		FILENAMEOUT = fmt.Sprintf("%s.%s.gz", INFILES[0], tag)
	case FILENAMEOUT == "" && BEDFILENAME != "":
		FILENAMEOUT = fmt.Sprintf("%s.%s.gz", BEDFILENAME, tag)
	case FILENAMEOUT == "":
		FILENAMEOUT = fmt.Sprintf("output.%s.gz", tag)
	}

	tStart := time.Now()

	switch {
	case COO || READINPEAK || CREATEBINMATRIX || TAIJI:
		switch {
		case CELLSIDFNAME == "" && !READINPEAK:
			log.Fatal("Error -xgi file must be provided!")
		case MERGEOUTPUTS:
			if len(INFILES) == 0 {
				log.Fatal("Error at least one input (-in) file must be provided!")
			}

			mergeMatFiles(INFILES)
		case BEDFILENAME == "":
			log.Fatal("Error at least one bed file must be provided!")
		case CREATEBINMATRIX:
			createBinSparseMatrix()
		case PEAKFILE == "" && !(COO || READINPEAK):
			log.Fatal("Error peak file -ygi (bed format) must be provided!")
		case READINPEAK:
			computeReadsInPeaksForCell()
		default:
			if SPLIT > 0 {
				createBoolSparseMatrixSplit()
			} else {
				createBoolSparseMatrix()
			}
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

	YGIDIM = utils.LoadPeaks(PEAKFILE)
	utils.CreatePeakIntervalTree()

	CELLIDCOUNT = make(map[string]int)

	fmt.Printf("init mutexes...\n")
	initMutexDict()
	fmt.Printf("init threading...\n")
	utils.InitIntervalDictsThreading(THREADNB)
	createReadInPeakOneFileThreading(BEDFILENAME)
	writeCellCounter(FILENAMEOUT)
}


func initIntSparseMatrix() {
	INTSPARSEMATRIX = make([]map[uint]int, XGIDIM)

	for _, pos := range CELLIDDICT {
		INTSPARSEMATRIX[pos] = make(map[uint]int)
	}
}


func createBoolSparseMatrix(){
	fmt.Printf("load indexes...\n")
	loadCellIDDict(CELLSIDFNAME)

	XGIDIM = len(CELLIDDICT)
	YGIDIM = utils.LoadPeaks(PEAKFILE)

	initIntSparseMatrix()
	utils.CreatePeakIntervalTree()
	launchBoolSparseMatrix(FILENAMEOUT, true)
}

func launchBoolSparseMatrix(filenameout string, writeTaijiHeader bool) {
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

	if TAIJI {
		writeIntMatrixToTaijiFile(filenameout, writeTaijiHeader)
	} else {
		writeIntMatrixToCOOFile(filenameout)
	}

}

func createBoolSparseMatrixSplit(){
	var count, nbsplit int
	var filenameout string
	var tmpfiles []string

	fmt.Printf("load indexes...\n")
	loadCellIDDict(CELLSIDFNAME)
	XGIDIM = len(CELLIDDICT)
	YGIDIM = utils.LoadPeaks(PEAKFILE)

	chunk := len(CELLIDDICT) / SPLIT
	celliddict := make([]string, len(CELLIDDICT))

	for cellID, pos := range CELLIDDICT {
		celliddict[pos] = cellID
	}

	CELLIDDICT = make(map[string]uint)

	utils.CreatePeakIntervalTree()

	for pos, cellID := range celliddict {
		CELLIDDICT[cellID] = uint(pos)

		count++

		if count > chunk {

			filenameout = fmt.Sprintf("%s.%d.tmp", FILENAMEOUT, nbsplit)
			fmt.Printf("#### Number of split: %d\n", nbsplit + 1)

			initIntSparseMatrix()
			launchBoolSparseMatrix(filenameout, nbsplit == 0)

			tmpfiles = append(tmpfiles, filenameout)

			count = 0
			nbsplit++
			CELLIDDICT = make(map[string]uint)
		}
	}

	//Finalizing last chunk
	filenameout = fmt.Sprintf("%s.%d.tmp", FILENAMEOUT, nbsplit)
	fmt.Printf("#### Number of split: %d\n", nbsplit + 1)
	initIntSparseMatrix()
	launchBoolSparseMatrix(filenameout, nbsplit == 0)
	tmpfiles = append(tmpfiles, filenameout)
	/////////////////////////////////////////////////////////////

	fmt.Printf("Concatenating tmp files...\n")
	cmd := fmt.Sprintf("cat %s > %s", strings.Join(tmpfiles, " "), FILENAMEOUT)
	utils.ExceCmd(cmd)

	fmt.Printf("file: %s created!\n", FILENAMEOUT)

	fmt.Printf("Removing tmp files...\n")
	cmd = fmt.Sprintf("rm %s", strings.Join(tmpfiles, " "))
	utils.ExceCmd(cmd)

}

func initMutexDict() {
	MUTEX = &sync.Mutex{}

	CELLMUTEXDICT = make(map[uint]*sync.Mutex)

	for _, pos := range CELLIDDICT {
		CELLMUTEXDICT[pos] = &sync.Mutex{}
	}
}


func writeIntMatrixToCOOFile(outfile string) {
	fmt.Printf("writing to output file...\n")

	var buffer bytes.Buffer
	var cellPos int
	var featPos uint
	var err error

	writer := utils.ReturnWriter(outfile)

	defer utils.CloseFile(writer)

	bufSize := 0

	for cellPos = range INTSPARSEMATRIX {
		for featPos = range INTSPARSEMATRIX[cellPos] {
			buffer.WriteString(strconv.Itoa(cellPos))
			buffer.WriteString(SEP)
			buffer.WriteString(strconv.Itoa(int(featPos)))
			buffer.WriteString(SEP)
			buffer.WriteString(strconv.Itoa(INTSPARSEMATRIX[cellPos][featPos]))
			buffer.WriteRune('\n')

			bufSize++

			if bufSize > 50000 {
				_, err = writer.Write(buffer.Bytes())
				utils.Check(err)
				buffer.Reset()
				bufSize = 0
			}
		}
	}

	_, err = writer.Write(buffer.Bytes())
	utils.Check(err)
	buffer.Reset()

	fmt.Printf("file: %s created!\n", outfile)
}

func writeCellCounter(outfile string) {
	var buffer bytes.Buffer
	var cellID string
	var err error

	writer := utils.ReturnWriter(outfile)

	defer utils.CloseFile(writer)

	bufSize := 0

	for cellID = range CELLIDCOUNT {
		if cellID == "" {
			continue
		}
		buffer.WriteString(cellID)
		buffer.WriteString(SEP)
		buffer.WriteString(strconv.Itoa(int(CELLIDCOUNT[cellID])))
		buffer.WriteRune('\n')

		bufSize++

		if bufSize > 50000 {
			_, err = writer.Write(buffer.Bytes())
			utils.Check(err)
			buffer.Reset()
			bufSize = 0

		}
	}

	_, err = writer.Write(buffer.Bytes())
			utils.Check(err)

	fmt.Printf("file: %s created!\n", outfile)
}

/*mergeMatFiles merge multiple COO output files*/
func mergeMatFiles(filenames []string) {
	fmt.Printf("creating xgi index..\n")
	loadCellIDDict(CELLSIDFNAME)

	if CREATEBINMATRIX {
		initBinSparseMatrix()
	} else {
		initIntSparseMatrix()
	}

	for _, filename := range filenames {
		fmt.Printf("merging file: %s\n", filename)

		if CREATEBINMATRIX {
			mergeCOOFloatMatFile(filename)
		} else {
			mergeCOOIntMatFile(filename)
		}
	}

	if TAIJI {
		if CREATEBINMATRIX {
			writeBinMatrixToTaijiFile(FILENAMEOUT)
		} else {
			writeIntMatrixToTaijiFile(FILENAMEOUT, true)
		}
	} else {
		if CREATEBINMATRIX {
			writeBinMatrixToCOOFile(FILENAMEOUT)
		} else {
			writeIntMatrixToCOOFile(FILENAMEOUT)
		}
	}
}


/*mergeCOOIntMatFile add one file to the matrix*/
func mergeCOOIntMatFile(filename string) {
	var scanner *bufio.Scanner
	var f *os.File
	var err error
	var split []string
	var xgi, ygi, value int

	scanner, f = utils.ReturnReader(filename, 0)

	defer utils.CloseFile(f)

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), SEP)
		xgi, err = strconv.Atoi(split[0])
		utils.Check(err)
		ygi, err = strconv.Atoi(split[1])
		utils.Check(err)
		value, err = strconv.Atoi(split[2])
		utils.Check(err)

		if ygi > YGIDIM {
			YGIDIM = ygi
		}

		switch USECOUNT {
		case true:
			INTSPARSEMATRIX[xgi][uint(ygi)] += value
		default:
			INTSPARSEMATRIX[xgi][uint(ygi)] = 1
		}

	}
}


/*mergeCOOFloatMatFile add one file to the matrix*/
func mergeCOOFloatMatFile(filename string) {
	var scanner *bufio.Scanner
	var f *os.File
	var err error
	var split []string
	var xgi, ygi int
	var value float64

	scanner, f = utils.ReturnReader(filename, 0)

	defer utils.CloseFile(f)

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), SEP)
		xgi, err = strconv.Atoi(split[0])
		utils.Check(err)
		ygi, err = strconv.Atoi(split[1])
		utils.Check(err)
		value, err = strconv.ParseFloat(split[2],64)
		utils.Check(err)

		if ygi > YGIDIM {
			YGIDIM = ygi
		}

		BINSPARSEMATRIX[uint(xgi)][uint(ygi)] += value

	}
}


/*createBoolSparseMatrixOneFileThreading ceate the bool Sparse Matrix for one bed file using multi-threading*/
func createBoolSparseMatrixOneFileThreading(bedfilename utils.Filename) {
	var nbReads uint
	var bufferLine1 [BUFFERSIZE]string
	var bufferLine2 [BUFFERSIZE]string
	var bufferPointer * [BUFFERSIZE]string

	isBuffer1 := true
	bufferPointer = &bufferLine1

	var bufferIt int
	var waiting sync.WaitGroup

	bedReader, file := bedfilename.ReturnReader(0)

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

func createReadInPeakOneFileThreading(bedfilename utils.Filename) {
	var nbReads uint
	var bufferLine1 [BUFFERSIZE]string
	var bufferLine2 [BUFFERSIZE]string
	var bufferPointer * [BUFFERSIZE]string

	var bufferIt int
	var waiting sync.WaitGroup

	isBuffer1 := true
	bufferPointer = &bufferLine1

	bedReader, file := bedfilename.ReturnReader(0)

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

		if len(intervals) == 0 {
			continue
		}

		CELLMUTEXDICT[cellPos].Lock()

		for _, int = range intervals {
			featPos = utils.PEAKIDDICT[utils.INTERVALMAPPING[int.ID()]]

			switch USECOUNT {
			case true:
				INTSPARSEMATRIX[cellPos][featPos]++
			default:
				INTSPARSEMATRIX[cellPos][featPos] = 1
			}
		}

		CELLMUTEXDICT[cellPos].Unlock()
	}
}

func createBoolSparseMatrixOneFile(bedfilename utils.Filename) {
	var line string
	var split []string
	var isInside bool
	var start, end int
	var err error
	var nbReads uint
	var intervals []interval.IntInterface
	var interval interval.IntInterface
	var cellPos,featPos uint

	bedReader, file := bedfilename.ReturnReader(0)

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

			switch USECOUNT {
			case true:
				INTSPARSEMATRIX[cellPos][featPos]++
			default:
				INTSPARSEMATRIX[cellPos][featPos] = 1
			}
		}
	}
}


func writeBinMatrixToTaijiFile(outfile string) {
	var buffer bytes.Buffer

	tStart := time.Now()

	sortedXgi := make([]string, len(BINSPARSEMATRIX))

	for xgi, pos := range CELLIDDICT {
		sortedXgi[pos] = xgi
	}

	writer := utils.ReturnWriter(FILENAMEOUT)
	defer utils.CloseFile(writer)

	buffer.WriteString("sparse matrix: ")
	buffer.WriteString(strconv.Itoa(len(sortedXgi)))
	buffer.WriteString(" x ")
	buffer.WriteString(strconv.Itoa(YGIDIM))
	buffer.WriteRune('\n')

	bufSize := 0

	var xgi int
	var ygi uint
	var value float64
	var err error

	for xgi = range BINSPARSEMATRIX {
		buffer.WriteString(sortedXgi[xgi])

		for ygi, value = range BINSPARSEMATRIX[xgi] {
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.Itoa(int(ygi)))
			buffer.WriteRune(',')

			if USECOUNT {
				buffer.WriteRune('1')

			} else
			{
				buffer.WriteString(strconv.FormatFloat(value, 'f', 8, 64))
			}

			bufSize++

			if bufSize > 10000 {
				_, err = writer.Write(buffer.Bytes())
				utils.Check(err)
				buffer.Reset()
				bufSize = 0
			}
		}

		buffer.WriteRune('\n')
	}

	_, err = writer.Write(buffer.Bytes())
	utils.Check(err)
	buffer.Reset()

	tDiff := time.Since(tStart)
	fmt.Printf("taiji formated matrix %s written in: %f s \n", FILENAMEOUT, tDiff.Seconds())

}


func writeIntMatrixToTaijiFile(outfile string, writeHeader bool) {
	var buffer bytes.Buffer

	tStart := time.Now()

	sortedXgi := make([]string, len(CELLIDDICT))

	for xgi, pos := range CELLIDDICT {
		sortedXgi[pos] = xgi
	}

	writer := utils.ReturnWriter(FILENAMEOUT)
	defer utils.CloseFile(writer)

	if writeHeader {
		buffer.WriteString("sparse matrix: ")
		buffer.WriteString(strconv.Itoa(len(sortedXgi)))
		buffer.WriteString(" x ")
		buffer.WriteString(strconv.Itoa(YGIDIM))
		buffer.WriteRune('\n')
	}

	bufSize := 0

	var ygi uint
	var xgi, value int
	var err error

	for xgi = range INTSPARSEMATRIX {
		buffer.WriteString(sortedXgi[xgi])

		for ygi, value = range INTSPARSEMATRIX[xgi] {
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.Itoa(int(ygi)))
			buffer.WriteRune(',')

			if USECOUNT {
				buffer.WriteRune('1')

			} else
			{
				buffer.WriteString(strconv.Itoa(value))
			}

			bufSize++

			if bufSize > 10000 {
				_, err = writer.Write(buffer.Bytes())
				utils.Check(err)
				buffer.Reset()
				bufSize = 0
			}
		}

		buffer.WriteRune('\n')
	}

	_, err = writer.Write(buffer.Bytes())
	utils.Check(err)
	buffer.Reset()

	tDiff := time.Since(tStart)
	fmt.Printf("taiji formated matrix %s written in: %f s \n", FILENAMEOUT, tDiff.Seconds())

}


/*loadCellIDDict load cell id to map[string] -> id <uint>*/
func loadCellIDDict(fname utils.Filename) {
	scanner, file := fname.ReturnReader(0)
	defer utils.CloseFile(file)
	var count uint
	var cellID string

	CELLIDDICT = make(map[string]uint)
	count = 0

	for scanner.Scan() {
		line := scanner.Text()
		cellID = strings.Split(line, "\t")[0]
		CELLIDDICT[cellID] = count
		count++
	}
}
