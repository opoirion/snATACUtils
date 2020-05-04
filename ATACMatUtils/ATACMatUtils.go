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
	"path"
	"sort"
	"math"
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

/*TOTALREADSCELL dict for total number of reads per cell map[cellID]total */
var TOTALREADSCELL []int

/*CELLIDREADINPEAK cell ID<->count reads in peaks  */
var CELLIDREADINPEAK map[string]int

/*CELLIDCOUNTALL cell ID<->count tot number of reads per cell */
var CELLIDCOUNTALL map[string]int

/*CELLIDDICT cell ID<->pos */
var CELLIDDICT map[string]uint

/*CELLIDDICTCOMP cell ID<->pos */
var CELLIDDICTCOMP []string

/*MUTEX global mutex */
var MUTEX *sync.Mutex

/*CELLMUTEXDICT feature pos<->sync.Mutex */
var CELLMUTEXDICT map[uint]*sync.Mutex

/*INTSPARSEMATRIX cell x feature sparse matrix  */
var INTSPARSEMATRIX []map[uint]int

/*FLOATSPARSEMATRIX cell x feature sparse matrix  */
var FLOATSPARSEMATRIX []map[uint]float64

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

/*YGISYMBOL use 4th columns  as ygi */
var YGISYMBOL bool

/*YGITOSYMBOL link ygi featurePos to symbol feature pos */
var YGITOSYMBOL []uint

/*YGISIZE size of YGI */
var YGISIZE []int

/*SYMBOLLIST ordered list of symbol */
var SYMBOLLIST []string

/*NORMTYPESTR string */
var NORMTYPESTR string

/*NORMTYPE normType */
var NORMTYPE normType

/*BUFFERSIZE buffer size for multithreading */
const BUFFERSIZE = 1000000

/*ALL bool indicating if all cells should be merged for count option */
var ALL bool

type normType string

const (
	fpkm normType = "fpkm"
	logfpkm normType = "logfpkm"
	rpm normType = "rpm"
	logrpm normType = "logrpm"
	simple normType = "simple"
	none normType = ""
)

func (t *normType) isValid() {
	switch *t {
	case none:
		return
	case rpm, simple, logrpm:
		NORM = true
		USECOUNT = true
	case fpkm, logfpkm:
		NORM = true
		USECOUNT = true

		if READINPEAK {
			panic(" Cannot use pfkm|logfpkm with -count ")
		}

	default:
		panic("Wrong -norm_type Normalisation type! possible value: rpm|fpkm|logrpm|logfpkm ")
	}
}

type mattype string

func (i *mattype) Set(mtype string) mattype {
	switch mtype {
	case "coo":
		*i = mattype(mtype)
	case "taiji":
		*i = mattype(mtype)
	default:
		panic(fmt.Sprintf("Matrix format unknown: %s", mtype))
	}

	return *i
}

func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `
#################### MODULE TO CREATE (cell x genomic region) SPARSE MATRIX ########################
"""Boolean / interger Peak matrix """
transform one (-bed) or multiple bed files into a sparse matrix
USAGE: ATACMatTools -coo -bed  <bedFile> -ygi <bedFile> -xgi <file> (-threads <int> -out <fname> -use_count -taiji -bed <bedFile2> -use_symbol)

"""Create a cell x bin matrix: -bin """
transform one (-bed) or multiple (use multiple -bed options) bed file into a bin (using float) sparse matrix. If ygi provided, reads intersecting these bin are ignored

USAGE: ATACMatTools -bin -bed  <bedFile> (optional -ygi <bedFile> -xgi <fname> -bin_size <int> -ygi_out <string> -norm -taiji -coo)

"""Count the number of reads in peaks for each cell: -count """
USAGE: ATACMatTools -count  -xgi <fname> -ygi <bedfile> -bed <bedFile> (optionnal: -out <fname> -norm -all)

"""Merge multiple matrices results into one output file: -merge """
It can be used to convert taiji to coo or coo to taiji formats.
USAGE: ATACMatTools -coo/taiji -merge -xgi <fname> -in <matrixFile1> -in <matrixFile2> ... (optional -bin -use_count -out <fname>)

USAGE for the -use_symbol option:
Optional use of the 4th column of -ygi as symbol.
This option is usefull to reduce the signal of multiple loci (i.e. multiple vectors) corresponding to the same function into one single vector per distinct symbol.
if used, the program will output the list of ordered symbol corresponding to the new ygi index.

`)
		 flag.PrintDefaults()
	}


	flag.IntVar(&SPLIT, "split", 0, "Split computation into n iterative chuncks (to reduce RAM usage for very large matrices)")
	flag.IntVar(&BINSIZE, "bin_size", 5000, "Size of the bin for bin matrix")
	flag.StringVar(&FILENAMEOUT, "out", "", "name of the output file")
	flag.Var(&BEDFILENAME, "bed", "name of the bed file")
	flag.Var(&INFILES, "in", "name of the input file(s)")
	flag.Var(&PEAKFILE, "ygi", "name of the bed file containing the region of interest( i.e. PEAK )")
	flag.BoolVar(&YGISYMBOL, "use_symbol", false, "Optional use of the 4th column of -ygi as symbol")
	flag.Var(&CELLSIDFNAME, "xgi", "name of the file containing the ordered list of cell IDs (one ID per line)")
	flag.StringVar(&SEP, "delimiter", "\t", "delimiter used to write the output file (default \t)")
	flag.StringVar(&YGIOUT, "ygi_out", "", "Write the output bin bed file in this specific file")
	flag.BoolVar(&TAIJI, "taiji", false,
		`Convert COO matrix file to sparse matrix format required for taiji-utils. See https://github.com/Taiji-pipeline/Taiji-utils/`)
	flag.BoolVar(&USECOUNT, "use_count", false,
		`Use read count instead of boolean value`)
	flag.BoolVar(&COO, "coo", false,
		`Use COO format as output`)
	flag.BoolVar(&ALL, "all", false,
		`Count the reads in peaks for the input bed file in its hole`)

	flag.BoolVar(&NORM, "norm", false, "Normalize raw count by dividing using the total number of reads per cell (equivalent to -norm_type simple)")
	flag.StringVar(&NORMTYPESTR, "norm_type", "", "Normalisation type to use: simple|rpm|log10(rpm)|fpkm. Simple divides the feature values by the total number of read count per cell (equivalent to norm)")

	flag.BoolVar(&CREATEBINMATRIX, "bin", false,
		`transform one (-bed) or multiple (use multiple -beds option) into a bin (using float) sparse matrix in COO format.`)
	flag.BoolVar(&MERGEOUTPUTS, "merge", false, `merge multiple matrices results into one output file`)
	flag.BoolVar(&READINPEAK, "count", false, `Count the number of reads in peaks for each cell`)
	flag.IntVar(&THREADNB, "threads", 1, "threads concurrency")
	flag.Parse()

	var tag string

	NORMTYPE = normType(NORMTYPESTR)

	if NORM && NORMTYPE == "" {
		NORMTYPE = "simple"
	}

	NORMTYPE.isValid()

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
				createIntSparseMatrixSplit()
			} else {
				createIntSparseMatrix()
			}
		}
	}

	tDiff := time.Since(tStart)
	fmt.Printf("done in time: %f s \n", tDiff.Seconds())
}


func loadSymbolFileWriteOutputSymbol() int {
	if !YGISYMBOL {
		return YGIDIM
	}

	if YGIDIM == 0 {
		panic(fmt.Sprintf("YGDIM is 0 when loading symbol"))
	}

	var symbol,line string
	var split []string

	SYMBOLLIST = []string{}

	scanner, file := PEAKFILE.ReturnReader(0)

	defer utils.CloseFile(file)

	symbolMap := make(map[int]string)
	symbolMapRev := make(map[string]int)
	symbolSet := make(map[string]bool)

	i := 0
	for scanner.Scan() {
		line = scanner.Text()
		split = strings.Split(line, "\t")

		if len(split) < 4 {
			panic(fmt.Sprintf(
				"Error when extracting symbol from %s cannot split line nb %d: %s in 4",
				PEAKFILE,
				i,
				line))
		}

		symbol = split[3]

		if !symbolSet[symbol] {
			symbolSet[symbol] = true
			SYMBOLLIST = append(SYMBOLLIST, symbol)
		}

		symbolMap[i] = symbol
		symbolMapRev[symbol] = i
		i++
	}

	sort.Strings(SYMBOLLIST)

	YGITOSYMBOL = make([]uint, YGIDIM)

	for indexNew, symbol  := range SYMBOLLIST {
		indexOld := symbolMapRev[symbol]
		YGITOSYMBOL[indexOld] = uint(indexNew)
	}

	fmt.Printf("symbol file loaded. New dim: %d\n", len(SYMBOLLIST))
	writeSymbol(SYMBOLLIST)

	return len(SYMBOLLIST)
}


func loadYgiSize() {

	switch NORMTYPE {
	case fpkm, logfpkm:
	default:
		return
	}

	var peak utils.Peak

	useSymbol := len(YGITOSYMBOL) != 0

	YGISIZE = make([]int, YGIDIM)

	if CREATEBINMATRIX {
		for pos := 0; pos < YGIDIM; pos ++ {
			YGISIZE[pos] = BINSIZE
		}

		return
	}

	for peakstr, pos := range utils.PEAKIDDICT {
		peak.StringToPeak(peakstr)

		if useSymbol {
			pos = YGITOSYMBOL[pos]
		}

		YGISIZE[pos] += peak.End - peak.Start
	}
}

func writeSymbol(SYMBOLLIST []string) {
	var buffer bytes.Buffer

	ext := path.Ext(FILENAMEOUT)
	filename := fmt.Sprintf("%s.symbol.ygi", FILENAMEOUT[:len(FILENAMEOUT) - len(ext)])

	writer := utils.ReturnWriter(filename)
	defer utils.CloseFile(writer)

	for _, symbol := range SYMBOLLIST {
		buffer.WriteString(symbol)
		buffer.WriteRune('\n')
	}

	writer.Write(buffer.Bytes())
	buffer.Reset()

	fmt.Printf("symbol file index written: %s\n", filename)
}

func computeReadsInPeaksForCell(){
	fmt.Printf("load indexes...\n")

	if CELLSIDFNAME != "" {
		loadCellIDDict(CELLSIDFNAME)
	}

	YGIDIM = utils.LoadPeaks(PEAKFILE)
	utils.CreatePeakIntervalTree()

	CELLIDREADINPEAK = make(map[string]int)

	if NORM {
		CELLIDCOUNTALL = make(map[string]int)
	}

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

	if NORM {
		TOTALREADSCELL = make([]int, XGIDIM)
	}
}


func initFloatSparseMatrix() {
	FLOATSPARSEMATRIX = make([]map[uint]float64, XGIDIM)
	BININDEX = make(map[binPos]uint)

	for _, pos := range CELLIDDICT {
		FLOATSPARSEMATRIX[pos] = make(map[uint]float64)
	}

	if NORM {
		TOTALREADSCELL = make([]int, XGIDIM)
	}
}


func createIntSparseMatrix(){
	fmt.Printf("load indexes...\n")
	loadCellIDDict(CELLSIDFNAME)

	XGIDIM = len(CELLIDDICT)
	YGIDIM = utils.LoadPeaks(PEAKFILE)
	YGIDIM = loadSymbolFileWriteOutputSymbol()

	initIntSparseMatrix()
	utils.CreatePeakIntervalTree()
	launchIntSparseMatrix(FILENAMEOUT, true)
}

func launchIntSparseMatrix(filenameout string, writeTaijiHeader bool) {
	switch{
	case THREADNB > 1:
		fmt.Printf("init mutexes...\n")
		initMutexDict()
		fmt.Printf("init threading...\n")
		utils.InitIntervalDictsThreading(THREADNB)
		fmt.Printf("launching sparse matrices creation...\n")
		createIntSparseMatrixOneFileThreading(BEDFILENAME)
	default:
		createIntSparseMatrixOneFile(BEDFILENAME)
	}

	if TAIJI {
		writeIntMatrixToTaijiFile(filenameout, writeTaijiHeader)
	} else {
		writeIntMatrixToCOOFile(filenameout)
	}

}

func createIntSparseMatrixSplit(){
	var count, nbsplit int
	var filenameout string
	var tmpfiles []string

	fmt.Printf("load indexes...\n")
	loadCellIDDict(CELLSIDFNAME)
	XGIDIM = len(CELLIDDICT)
	YGIDIM = utils.LoadPeaks(PEAKFILE)
	YGIDIM = loadSymbolFileWriteOutputSymbol()

	chunk := XGIDIM / SPLIT
	celliddict := make([]string, XGIDIM)

	for cellID, pos := range CELLIDDICT {
		celliddict[pos] = cellID
	}

	CELLIDDICT = make(map[string]uint)

	utils.CreatePeakIntervalTree()
	ext := path.Ext(FILENAMEOUT)

	for pos, cellID := range celliddict {
		CELLIDDICT[cellID] = uint(pos)

		count++

		if count > chunk {

			filenameout = fmt.Sprintf("%s.%d.tmp%s", FILENAMEOUT[:len(FILENAMEOUT) - len(ext)], nbsplit, ext)
			fmt.Printf("#### Number of split: %d\n", nbsplit + 1)

			initIntSparseMatrix()
			launchIntSparseMatrix(filenameout, nbsplit == 0)

			tmpfiles = append(tmpfiles, filenameout)

			count = 0
			nbsplit++
			CELLIDDICT = make(map[string]uint)
		}
	}

	//Finalizing last chunk
	filenameout = fmt.Sprintf("%s.%d.tmp%s", FILENAMEOUT[:len(FILENAMEOUT) - len(ext)], nbsplit, ext)
	fmt.Printf("#### Number of split: %d\n", nbsplit + 1)
	initIntSparseMatrix()
	launchIntSparseMatrix(filenameout, nbsplit == 0)
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
	loadYgiSize()

	var buffer bytes.Buffer
	var cellPos int
	var featPos uint
	var err error
	var normedValue float64

	writer := utils.ReturnWriter(outfile)

	defer utils.CloseFile(writer)

	bufSize := 0

	for cellPos = range INTSPARSEMATRIX {
		for featPos = range INTSPARSEMATRIX[cellPos] {

			buffer.WriteString(strconv.Itoa(cellPos))
			buffer.WriteString(SEP)
			buffer.WriteString(strconv.Itoa(int(featPos)))
			buffer.WriteString(SEP)

			if NORM {
				normedValue = normValue(
					INTSPARSEMATRIX[cellPos][featPos],
					cellPos, int(featPos))
				buffer.WriteString(strconv.FormatFloat(normedValue, 'f', 7, 64))
			} else {
				buffer.WriteString(strconv.Itoa(INTSPARSEMATRIX[cellPos][featPos]))
			}

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

func writeFloatMatrixToCOOFile(outfile string) {
	fmt.Printf("writing to output file...\n")

	var buffer bytes.Buffer
	var cellPos int
	var featPos uint
	var err error
	var normedValue float64

	writer := utils.ReturnWriter(outfile)

	defer utils.CloseFile(writer)

	bufSize := 0

	for cellPos = range FLOATSPARSEMATRIX {
		for featPos = range FLOATSPARSEMATRIX[cellPos] {

			buffer.WriteString(strconv.Itoa(cellPos))
			buffer.WriteString(SEP)
			buffer.WriteString(strconv.Itoa(int(featPos)))
			buffer.WriteString(SEP)

			normedValue = FLOATSPARSEMATRIX[cellPos][featPos]
			buffer.WriteString(strconv.FormatFloat(normedValue, 'f', 7, 64))

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

func normValue(value int, cellID, featID int) (valueFloat float64) {
	valueFloat  = float64(value)

	switch NORMTYPE {
	case simple:
		valueFloat = valueFloat / float64(TOTALREADSCELL[cellID])
	case rpm:
		valueFloat = valueFloat * 1e6 / float64(TOTALREADSCELL[cellID])
	case logrpm:
		valueFloat = math.Log((1 + valueFloat) * 1e6 /
			float64(TOTALREADSCELL[cellID]))
	case fpkm:
		valueFloat = valueFloat * 1e9 /
			(float64(TOTALREADSCELL[cellID]) * float64(YGISIZE[featID]))
	case logfpkm:
		valueFloat = math.Log((1 + valueFloat) * 1e9 /
			(float64(TOTALREADSCELL[cellID]) * float64(YGISIZE[featID])))
	}

	return valueFloat
}


func normValueForReadInPeaks(value float64, cellIDstr string) float64 {
	switch NORMTYPE {
	case simple:
		value = float64(CELLIDREADINPEAK[cellIDstr]) / float64(CELLIDCOUNTALL[cellIDstr])
	case rpm:
		value = float64(CELLIDREADINPEAK[cellIDstr]) * 1e6 / float64(CELLIDCOUNTALL[cellIDstr])
	case logrpm:
		value = float64(CELLIDREADINPEAK[cellIDstr]) * 1e6 / float64(CELLIDCOUNTALL[cellIDstr])
		value = math.Log10(value)
	}

	return value
}

func writeCellCounter(outfile string) {
	var buffer bytes.Buffer
	var cellID string
	var err error
	var value int
	var valueFloat float64

	writer := utils.ReturnWriter(outfile)

	defer utils.CloseFile(writer)

	bufSize := 0


	for cellID, value = range CELLIDREADINPEAK {
		if cellID == "" {
			continue
		}
		buffer.WriteString(cellID)
		buffer.WriteString(SEP)

		if NORM {
			valueFloat = normValueForReadInPeaks(float64(value), cellID)

			buffer.WriteString(strconv.FormatFloat(valueFloat, 'f', 7, 64))
		} else {
			buffer.WriteString(strconv.Itoa(value))
		}


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
	XGIDIM = len(CELLIDDICT)

	if CREATEBINMATRIX {
		initFloatSparseMatrix()
	} else {
		initIntSparseMatrix()
	}

	for _, filename := range filenames {
		mtype := findMatrixFormat(filename)
		fmt.Printf("merging file: matrix %s with format: %s \n", filename, mtype)

		if CREATEBINMATRIX {
			mergeFloatMatFile(filename, mtype)
		} else {
			mergeIntMatFile(filename, mtype)
		}
	}

	if TAIJI {
		if CREATEBINMATRIX {
			writeFloatMatrixToTaijiFile(FILENAMEOUT)
		} else {
			writeIntMatrixToTaijiFile(FILENAMEOUT, true)
		}
	} else {
		if CREATEBINMATRIX {
			writeFloatMatrixToCOOFile(FILENAMEOUT)
		} else {
			writeIntMatrixToCOOFile(FILENAMEOUT)
		}
	}
}

/*mergeIntMatFile add one file to the matrix*/
func mergeIntMatFile(filename string, mtype mattype) {
	switch mtype {
	case "coo":
		mergeIntMatFileFromCOO(filename)
	case "taiji":
		mergeIntMatFileFromTaiji(filename)
	}
}

/*mergeFloatMatFile add one file to the matrix*/
func mergeFloatMatFile(filename string, mtype mattype) {
	switch mtype {
	case "coo":
		mergeFloatMatFileFromCOO(filename)
	case "taiji":
		mergeFloatMatFileFromTaiji(filename)
	}
}

/*findMatrixFormat find matrix format */
func findMatrixFormat(filename string) mattype {
	scanner, f := utils.ReturnReader(filename, 0)

	scanner.Scan()
	firstLine := scanner.Text()
	utils.CloseFile(f)

	var err1, err2, err3 error

	split := strings.Split(firstLine, "\t")

	if len(split) == 3 {
		_, err1 = strconv.Atoi(split[0])
		_, err2 = strconv.Atoi(split[1])
		_, err3 = strconv.ParseFloat(split[2], 64)

		if err1 == nil && err2 == nil && err3 == nil {
			return mattype("coo")
		}
	}

	if strings.Contains(firstLine, "parse matrix:") && strings.Contains(firstLine, " x ") {
		return mattype("taiji")
	}

	panic(fmt.Sprintf("Matrix header:% s does not match any matrix type", firstLine))

	return mattype("unknown")
}

/*mergeIntMatFileFromCOO add one file to the matrix*/
func mergeIntMatFileFromCOO(filename string) {
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


func getTaijiMatDim(filename string) int {
	scanner, f := utils.ReturnReader(filename, 0)
	scanner.Scan()

	firstLine := scanner.Text()
	utils.CloseFile(f)

	ygidim, err := strconv.Atoi(strings.Split(firstLine, " x ")[1])
	utils.Check(err)

	return ygidim

}

/*mergeIntMatFileFromCOO add one file to the matrix*/
func mergeIntMatFileFromTaiji(filename string) {
	var scanner *bufio.Scanner
	var f *os.File
	var err error
	var split, ygiSplit []string
	var ygiUnit string
	var xgi uint
	var ygi, ygidim int
	var value float64
	var isInside bool

	ygidim = getTaijiMatDim(filename)
	scanner, f = utils.ReturnReader(filename, 0)
	defer utils.CloseFile(f)

	scanner.Buffer([]byte{}, ygidim / 3 * ygidim / 3)
	scanner.Scan()

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), "\t")

		if xgi, isInside = CELLIDDICT[split[0]]; !isInside {
			fmt.Printf("Barcode %s is not in %s", split[0], CELLSIDFNAME)
			continue
		}

		for _, ygiUnit = range split[1:] {
			ygiSplit = strings.Split(ygiUnit, ",")
			ygi, err = strconv.Atoi(ygiSplit[0])
			utils.Check(err)
			value, err = strconv.ParseFloat(ygiSplit[1], 64)
			utils.Check(err)

			if ygi > YGIDIM {
				YGIDIM = ygi
			}

			switch USECOUNT {
			case true:
				INTSPARSEMATRIX[xgi][uint(ygi)] += int(value)
			default:
				INTSPARSEMATRIX[xgi][uint(ygi)] = 1
			}

		}
	}
}


/*mergeFLloattMatFileFromCOO add one file to the matrix*/
func mergeFloatMatFileFromTaiji(filename string) {
	var scanner *bufio.Scanner
	var f *os.File
	var err error
	var split, ygiSplit []string
	var ygiUnit string
	var xgi uint
	var ygi, ygidim int
	var value float64
	var isInside bool

	ygidim = getTaijiMatDim(filename)
	scanner, f = utils.ReturnReader(filename, 0)
	defer utils.CloseFile(f)

	scanner.Buffer([]byte{}, ygidim / 3 * ygidim / 3)
	scanner.Scan()

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), "\t")

		if xgi, isInside = CELLIDDICT[split[0]]; !isInside {
			fmt.Printf("Barcode %s is not in %s", split[0], CELLSIDFNAME)
			continue
		}

		for _, ygiUnit = range split[1:] {
			ygiSplit = strings.Split(ygiUnit, ",")
			ygi, err = strconv.Atoi(ygiSplit[0])
			utils.Check(err)
			value, err = strconv.ParseFloat(ygiSplit[1], 64)
			utils.Check(err)

			if ygi > YGIDIM {
				YGIDIM = ygi
			}

			FLOATSPARSEMATRIX[xgi][uint(ygi)] += value
		}
	}
}


/*mergeFloatMatFileFromCOO add one file to the matrix*/
func mergeFloatMatFileFromCOO(filename string) {
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
		value, err = strconv.ParseFloat(split[2], 64)
		utils.Check(err)

		if ygi > YGIDIM {
			YGIDIM = ygi
		}

		FLOATSPARSEMATRIX[uint(xgi)][uint(ygi)] += value

	}
}


/*createIntSparseMatrixOneFileThreading ceate the bool Sparse Matrix for one bed file using multi-threading*/
func createIntSparseMatrixOneFileThreading(bedfilename utils.Filename) {
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
				go updateIntSparseMatrixOneThread(
					bufferPointer , bufferStart, bufferStop, i, &waiting)

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
		updateIntSparseMatrixOneThread(bufferPointer , 0, bufferIt, 0, &waiting)
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
				go updateReadInPeakThread(bufferPointer,
					bufferStart, bufferStop, i, &waiting)

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
	var start, end, count int
	var err error
	var cellcount []string
	var cellIDstr string

	isCellsID := true

	if CELLSIDFNAME == "" {
		isCellsID = false
	}

	if NORM {
		cellcount = make([]string, bufferStop - bufferStart)
	}

	var intervals []interval.IntInterface

	for i := bufferStart; i < bufferStop;i++ {

		split = strings.Split(bufferLine[i], "\t")
		if ALL {
			cellIDstr = "all"
		} else {
			cellIDstr = split[3]
		}

		if isCellsID {
			if _, isInside = CELLIDDICT[split[3]];!isInside {
				continue
			}
		}

		if NORM {
			cellcount[count] = cellIDstr
			count++
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

		if len(intervals) > 0 {
			CELLIDREADINPEAK[cellIDstr]++
		}

		if NORM {
			for _,cellIDstr = range cellcount[:count] {
				CELLIDCOUNTALL[cellIDstr]++
			}

			count = 0
		}

		MUTEX.Unlock()
	}

	MUTEX.Lock()

	//process last reads for norm
	if NORM {
		for _,cellIDstr = range cellcount[:count] {
			CELLIDCOUNTALL[cellIDstr]++
		}
	}

	MUTEX.Unlock()
}

func updateIntSparseMatrixOneThread(
	bufferLine * [BUFFERSIZE]string,
	bufferStart,
	bufferStop,
	threadnb int,
	waiting * sync.WaitGroup) {
	defer waiting.Done()
	var split []string
	var isInside bool
	var start, end int
	var err error

	var intervals []interval.IntInterface
	var inter interval.IntInterface
	var cellPos,featPos, cellPosRead uint
	var totalreadscell []uint

	if NORM {
		totalreadscell = make([]uint, BUFFERSIZE)
	}

	useSymbol := len(YGITOSYMBOL) != 0

	for i := bufferStart; i < bufferStop;i++ {

		split = strings.Split(bufferLine[i], "\t")

		if cellPos, isInside = CELLIDDICT[split[3]];!isInside {
			continue
		}

		if NORM {
			totalreadscell[i] = cellPos + 1
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

		for _, inter = range intervals {
			featPos = utils.PEAKIDDICT[utils.INTERVALMAPPING[inter.ID()]]

			if useSymbol {
				featPos = YGITOSYMBOL[featPos]
			}

			switch USECOUNT {
			case true:
				INTSPARSEMATRIX[cellPos][featPos]++
			default:
				INTSPARSEMATRIX[cellPos][featPos] = 1
			}
		}

		CELLMUTEXDICT[cellPos].Unlock()
	}

	if NORM {
		BININDEXMUTEX.Lock()
		for _, cellPosRead = range totalreadscell[bufferStart:bufferStop] {

			if cellPosRead > 0 {
				TOTALREADSCELL[cellPosRead-1]++
			}
		}
		BININDEXMUTEX.Unlock()
	}
}

func createIntSparseMatrixOneFile(bedfilename utils.Filename) {
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

		if NORM {
			TOTALREADSCELL[cellPos]++
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


func writeFloatMatrixToTaijiFile(outfile string) {
	var buffer bytes.Buffer

	tStart := time.Now()

	sortedXgi := make([]string, XGIDIM)

	for xgi, pos := range CELLIDDICT {
		sortedXgi[pos] = xgi
	}

	writer := utils.ReturnWriter(FILENAMEOUT)
	defer utils.CloseFile(writer)

	buffer.WriteString("Sparse matrix: ")
	buffer.WriteString(strconv.Itoa(XGIDIM))
	buffer.WriteString(" x ")
	buffer.WriteString(strconv.Itoa(YGIDIM))
	buffer.WriteRune('\n')

	bufSize := 0

	var xgi int
	var ygi uint
	var value float64
	var err error

	for xgi = range FLOATSPARSEMATRIX {
		buffer.WriteString(sortedXgi[xgi])

		for ygi, value = range FLOATSPARSEMATRIX[xgi] {
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.Itoa(int(ygi)))
			buffer.WriteRune(',')

			if NORM {
				buffer.WriteString(strconv.FormatFloat(value, 'f', 7, 64))
			} else
			{
				buffer.WriteRune('1')
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


func writeIntMatrixToTaijiFile(filenameout string, writeHeader bool) {
	var buffer bytes.Buffer
	var floatValue float64

	tStart := time.Now()

	sortedXgi := make([]string, XGIDIM)

	for xgi, pos := range CELLIDDICT {
		sortedXgi[pos] = xgi
	}

	writer := utils.ReturnWriter(filenameout)
	defer utils.CloseFile(writer)
	loadYgiSize()

	if writeHeader {
		buffer.WriteString("Sparse matrix: ")
		buffer.WriteString(strconv.Itoa(XGIDIM))
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

			if NORM {
				floatValue = normValue(
					value, xgi, int(ygi))
				buffer.WriteString(strconv.FormatFloat(
					floatValue, 'f', 7, 64))
			} else {
				buffer.WriteRune('1')
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
	fmt.Printf("taiji formated matrix %s written in: %f s \n", filenameout, tDiff.Seconds())

}


/*loadCellIDDict load cell id to map[string] -> id <uint>*/
func loadCellIDDict(fname utils.Filename) {
	scanner, file := fname.ReturnReader(0)
	defer utils.CloseFile(file)
	var count uint
	var cellID string
	var isInside bool

	CELLIDDICT = make(map[string]uint)
	count = 0

	for scanner.Scan() {
		line := scanner.Text()
		cellID = strings.Split(line, "\t")[0]

		if _, isInside = CELLIDDICT[cellID];isInside {
			panic(fmt.Sprintf("At line: %d cellID: %s is present twice in: %s",
				count, cellID, fname))
		}

		CELLIDDICT[cellID] = count
		count++
	}

	CELLIDDICTCOMP = make([]string, len(CELLIDDICT))

	for cellID, count = range CELLIDDICT {
		CELLIDDICTCOMP[count] = cellID
	}
}
