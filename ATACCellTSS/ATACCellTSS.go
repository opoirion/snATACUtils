/* Create excecutable to compute individual TSS for each cell */

package main


import(
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"flag"
	"os"
	"io"
	"fmt"
	"log"
	"sync"
	"strings"
	"strconv"
	"github.com/biogo/store/interval"
	"time"
	"bytes"
	"path"
	"math"
)


//////////////// INSTANCES ////////////////////////////////////

/*BASECOVERAGE map[cell][tss relative pos]count */
var BASECOVERAGE [][]int

/*FLANKCOVERAGE map[cell]count */
var FLANKCOVERAGE []int

/*CELLTSS map[cell]tss */
var CELLTSS map[string]float64

/*CELLCLUSTERDICT map[cell]clusterID*/
var CELLCLUSTERDICT map[string]int

/*CELLDICT map[cell]int */
var CELLDICT map[string]int

/*BUFFERSIZE buffer size for multithreading */
const BUFFERSIZE = 100000

/*MUTEX global mutex*/
var MUTEX sync.Mutex

/*BUFFERARRAY buffer array storing lines of bed files*/
var BUFFERARRAY[][BUFFERSIZE]string

/*FLANKCOVERAGEMAT mat of array for flank cov when creating matrix
genomic region -> clID (int) -> flank coverage vector
*/
var FLANKCOVERAGEMAT map[uintptr][]int

/*BASECOVERAGEMAT mat of array for base cov when creating matrix
genomic region -> clID (int) -> base coverage vector
 */
var BASECOVERAGEMAT map[uintptr][][]int

/*BUFFERLENGTH effective size of buffer int*/
var BUFFERLENGTH int

/*ORIENTATIONDICT mapping peaks to orientation (e.g. +/1) */
var ORIENTATIONDICT map[uint]string
////////////////////////////////////////////////////////////////


//////////////////// INPUT VARIABLES ////////////////////////////////////////

/*MATRIXBINSIZE bin size of the matrix */
var MATRIXBINSIZE uint

/*BEDFILENAME bed file name (input) */
var BEDFILENAME utils.Filename

/*PEAKFILE bed file name (input) */
var PEAKFILE utils.Filename

/*TSSFILE bed file name (input) */
var TSSFILE utils.Filename

/*CELLSIDFNAME file name file with ordered cell IDs (one ID per line) */
var CELLSIDFNAME utils.Filename

/*CLUSTERFNAME file name file with cluster<TAB>cell IDs for each line */
var CLUSTERFNAME utils.Filename

/*FILENAMEOUT  output file name output */
var FILENAMEOUT string

/*FLANKSIZE int both at the end and begining of TSS region*/
var FLANKSIZE int

/*TSSREGION int both at the end and begining of TSS region*/
var TSSREGION int

/*SHIFTREAD shift read in the 3'->5' orientation*/
var SHIFTREAD int

/*SMOOTHINGWINDOW int smoothing window size*/
var SMOOTHINGWINDOW int

/*TSSFLANKSEARCH tss flank size*/
var TSSFLANKSEARCH int

/*THREADNB number of threads for reading the bam file */
var THREADNB int

/*USEMIDDLE bool*/
var USEMIDDLE bool

/*COLSEQID int*/
var COLSEQID int

/*REFCOLSEQID int*/
var REFCOLSEQID int

/*ALL bool*/
var ALL bool

/*STDOUT write to stdout*/
var STDOUT bool

/*CREATEMAT create TSS enrichment matrix for plotting (using plotHeatmap from deepTools) */
var CREATEMAT bool
////////////////////////////////////////////////////////////////////


func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `
USAGE: ATACCellTSS
                       -bed <filename>
                       -ygi/tss <filename>
                       -xgi <filename>
##### optional ####
                       -out <string>
                       -flank <int>
                       -smoothing <int>
                       -boundary <int>
                       -cluster <filename>
                       -flank_size
                       -threads <int>
                       -col_seqID <int>
                       -col_refID <int>
                       -xgi <filename>
                       -create_TSS_matrix
                       -bin_size <int>

if -cluster is provided, TSS is computed per cluster and -xgi argument is ignored. THe cluster file should contain cluster and cell ID with the following structure for each line: clusterID<TAB>cellID\n

if -create_TSS_matrix is provided, the program will output in addition a matrix file containing a TSS enrichment matrix file for each group (if -cluster is provided) or a global matrix (if -all is provided)

`)
		 flag.PrintDefaults()
	}

	flag.Var(&BEDFILENAME, "bed", "name of the bed file")
	flag.Var(&PEAKFILE, "ygi", "name of the bed file containing the TSS regions (assuming TSS is at the center)")
	flag.Var(&TSSFILE, "tss", "name of the bed file containing the TSS position (alternative to ygi)")
	flag.Var(&CELLSIDFNAME, "xgi", "name of the file containing the ordered list of cell IDs (one ID per line)")
	flag.Var(&CLUSTERFNAME, "cluster", "name of the file containing the cluster<->cellID (<TAB> separated). If cluster is provided, then -xgi is ignored")
	flag.StringVar(&FILENAMEOUT, "out", "", "name of the output file")
	flag.IntVar(&FLANKSIZE, "flank", 100, "flank size at the end and begining of the TSS regions")
	flag.IntVar(&COLSEQID, "col_seqID", -1, "If > 0, use this column as additional sequence ID (e.g. orientation +/-)")
	flag.IntVar(&REFCOLSEQID, "col_refID", -1, "If > 0, use this column as additional sequence ID (e.g. orientation +/-) for reference features")
	flag.IntVar(&TSSREGION, "boundary", 2000, "TSS boundary size at the end and begining of the TSS (used only when -tss is provided)")
	flag.IntVar(&SMOOTHINGWINDOW, "smoothing", 50, "Smoothing window size")
	flag.IntVar(&TSSFLANKSEARCH, "tss_flank", 50, "search hightest TSS values to define TSS score using this flank size arround TSS")
	flag.IntVar(&THREADNB, "threads", 1, "threads concurrency")
	flag.IntVar(&SHIFTREAD, "shift_reads", 0, "shift bed reads by adding X bp. If a \"-\" orientation is given, this number is substracted ")
	flag.UintVar(&MATRIXBINSIZE, "bin_size", 10, "Bin size to average scores when constructing TSS matrix")
	flag.BoolVar(&USEMIDDLE, "use_middle", false, "Use the middle of the peak to determine the TSS ")
	flag.BoolVar(&STDOUT, "stdout", false, `write to stdout`)
	flag.BoolVar(&CREATEMAT, "create_TSS_matrix", false, `create TSS enrichment matrix for plotting (using plotHeatmap from deepTools)`)
	flag.BoolVar(&ALL, "all", false, "Compute the general TSS ")
	flag.Parse()

	if MATRIXBINSIZE == 0 {
		panic("Error -bin_size should be higher than 0")
	}

	switch {
	case BEDFILENAME == "":
		log.Fatal("-bed must be provided!")
	case PEAKFILE == "" && TSSFILE == "":
		log.Fatal("either -ygi or -tss must be provided!")

	case FILENAMEOUT == "":
		FILENAMEOUT = findFileNameOut()
	}

	if CREATEMAT && !ALL && CLUSTERFNAME == "" {
		fmt.Printf("Error with -create_TSS_matrix option. Either -all or -cluster shpuld be provided\n")
	}

	MUTEX = sync.Mutex{}

	if PEAKFILE != "" {
		if USEMIDDLE {
			loadTSS(PEAKFILE)
		} else {
			_, ORIENTATIONDICT = utils.LoadPeaksAndTrimandReturnOrientation(
				PEAKFILE, REFCOLSEQID)
		}

	} else {
		loadTSS(TSSFILE)
	}

	switch {

	case CLUSTERFNAME != "":
		loadClusterFile()
	case ALL && CELLSIDFNAME == "":
		//Nothing
	case CELLSIDFNAME != "":
		CELLDICT = utils.LoadCellDictsToIndex(CELLSIDFNAME)
	default:
		CELLDICT = utils.LoadCellDictsFromBedFileToIndex(BEDFILENAME)
	}

	if ALL {
		makeClusterFileForAll()
	}

	utils.CreatePeakIntervalTree()
	utils.InitIntervalDictsThreading(THREADNB)

	initDicts()

	scanBedFile()
	normaliseCountMatrix()
	writeCellTSSScore()

	if CREATEMAT {
		writeTSSMatrix()
	}
}


func findFileNameOut() (foutname string) {

	switch {
	case CLUSTERFNAME != "":
		foutname = CLUSTERFNAME.String()
	case CELLSIDFNAME == "":
		foutname = BEDFILENAME.String()
	default:
		foutname = CELLSIDFNAME.String()
	}

	ext := path.Ext(foutname)

	var trueExt string

	switch {
	case ALL:
		trueExt = ".tss"
	case CLUSTERFNAME != "":
		trueExt = ".tss_per_cluster"
	default:
		trueExt = ".tss_per_cell"
	}

	foutname = fmt.Sprintf("%s%s",
		foutname[:len(foutname)-len(ext)], trueExt)

	return foutname
}

func loadClusterFile() {
	var split []string
	var clusterID, index int
	var cell, cluster string
	var isInside bool

	CELLDICT = make(map[string]int)
	CELLCLUSTERDICT = make(map[string]int)

	scanner, file := CLUSTERFNAME.ReturnReader(0)
	defer utils.CloseFile(file)

	scanner.Scan()

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), "\t")
		cell, cluster = split[0], split[1]

		if clusterID, isInside = CELLDICT[cluster]; !isInside {
			CELLDICT[cluster] = index
			clusterID = index
			index++
		}

		CELLCLUSTERDICT[cell] = clusterID
	}

}


func makeClusterFileForAll() {
	var cell string

	CELLCLUSTERDICT = make(map[string]int)

	for cell = range CELLDICT {

		CELLCLUSTERDICT[cell] = 0
	}

	CELLDICT = map[string]int{"all":0}
}


func loadTSS(tssFile utils.Filename) {

	scanner, file := tssFile.ReturnReader(0)
	defer utils.CloseFile(file)

	var count uint
	var start int
	var end int
	var buffer bytes.Buffer
	var split []string
	var err error
	var chro string

	utils.PEAKIDDICT = make(map[string]uint)
	count = 0

	for scanner.Scan() {
		split = strings.Split(
			strings.Trim(scanner.Text(), "\n"), "\t")
		chro = strings.TrimPrefix(split[0], "chr")

		buffer.WriteString(chro)
		buffer.WriteRune('\t')

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		if USEMIDDLE {
			end, err = strconv.Atoi(split[2])
			utils.Check(err)

			start = (start + end) / 2
		}

		buffer.WriteString(strconv.Itoa(start - TSSREGION))
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.Itoa(start + TSSREGION))

		if len(split) >= 4 {
			buffer.WriteRune('\t')
			buffer.WriteString(split[3])
		}

		buffer.WriteRune('\n')

		utils.PEAKIDDICT[buffer.String()] = count
		buffer.Reset()
		count++
	}
}

func initCoverageMats(intervalID uintptr) {

	if _, isInside := FLANKCOVERAGEMAT[intervalID];isInside {
		return
	}

	var index int

	FLANKCOVERAGEMAT[intervalID] = make([]int, len(CELLDICT))
	BASECOVERAGEMAT[intervalID] = make([][]int, len(CELLDICT))

	for _, index = range CELLDICT {
		BASECOVERAGEMAT[intervalID][index] = make([]int, 2 * TSSREGION)
	}

}

func initDicts() {
	var index int

	FLANKCOVERAGE = make([]int, len(CELLDICT))
	BASECOVERAGE = make([][]int, len(CELLDICT))
	CELLTSS = make(map[string]float64)

	if CREATEMAT {
		FLANKCOVERAGEMAT = make(map[uintptr][]int, len(CELLDICT))
		BASECOVERAGEMAT = make(map[uintptr][][]int, len(CELLDICT))

	}

	TSSREGION = assertPeakRegionHaveSameLengths()

	// if PEAKFILE != "" {
	// 	TSSREGION = assertPeakRegionHaveSameLengths()
	// }

	BUFFERLENGTH = 2 * TSSREGION - 2 * FLANKSIZE

	if 2 * FLANKSIZE >= TSSREGION {
		log.Fatal(fmt.Sprintf("Error 2 * %d (flank size) greater than %d (TSS region) \n",
			FLANKSIZE, TSSREGION))
	}

	if TSSREGION <= TSSFLANKSEARCH {
		log.Fatal(fmt.Sprintf("Error TSS region (-boundary):%d must be greater than tss flank search window (-tss_flank):%d \n",
			TSSREGION, TSSFLANKSEARCH))
	}

	fmt.Printf("TSS region length: %d actual screening window size for TSS score: %d\n",
		TSSREGION, BUFFERLENGTH)

	for _, index = range CELLDICT {
		BASECOVERAGE[index] = make([]int, BUFFERLENGTH)
	}
}


func assertPeakRegionHaveSameLengths() (tssregion int){
	var peak string
	var split []string
	var start, end int
	var err error

	tssregion = 0

	for peak = range utils.PEAKIDDICT {
		split = strings.Split(peak, "\t")

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(strings.Trim(split[2], "\n"))
		utils.Check(err)

		if tssregion == 0 {
			tssregion = end - start
		} else {
			if end - start != tssregion {
				log.Fatal(fmt.Sprintf(
					`Constant TSS tssregion error!
Error for TSS region: %s, expected to have a size of %d, found %d instead \n
To use peaks with different sizes, please use the -use_middle option
`,
					peak, tssregion, end - start))
			}
		}
	}

	return tssregion / 2
}


func scanBedFile() {
	var count, thread int
	var waiting sync.WaitGroup
	var freeThreads = make(chan int, THREADNB)

	for i := 0; i < THREADNB; i++ {
		freeThreads <- i
	}

	thread =<-freeThreads

	BUFFERARRAY = make([][BUFFERSIZE]string, THREADNB)

	scanner, file := BEDFILENAME.ReturnReader(0)
	defer utils.CloseFile(file)

	tStart := time.Now()

	for scanner.Scan() {
		BUFFERARRAY[thread][count] = scanner.Text()
		count++

		if count >= BUFFERSIZE {
			waiting.Add(1)
			go processOneBuffer(&BUFFERARRAY[thread], thread, count, freeThreads, &waiting)
			count = 0
			thread = <- freeThreads
		}
	}

	waiting.Add(1)
	go processOneBuffer(&BUFFERARRAY[thread], thread, count, freeThreads, &waiting)
	waiting.Wait()

	tDiff := time.Since(tStart)
	fmt.Printf("Scanning bed file done in time: %f s \n", tDiff.Seconds())
}


func processOneBuffer(bufferarray *[BUFFERSIZE]string,
	thread, limit int,
	freeThreads chan int, waiting *sync.WaitGroup) {

	var split []string
	var chro string
	var start, end, j, index, cellID int
	var err error
	var isInside bool
	var intervals []interval.IntInterface
	var inter interval.IntInterface
	var itrg interval.IntRange
	var celldict map[string]int
	var cellIDis0 bool
	var intervalID uintptr

	var basecoverage [][]int
	var flankcoverage []int

	var readOrientation, refOrientation string
	useOrientation := COLSEQID > -1
	doShiftting := SHIFTREAD != 0

	defer waiting.Done()

	// Check if TSS is computed per cell or per cluster
	switch {
	case ALL:
		celldict = CELLCLUSTERDICT
	case CLUSTERFNAME != "":
		celldict = CELLCLUSTERDICT
	default:
		celldict = CELLDICT
	}

	indexLimit := 2 * TSSREGION - 2 * FLANKSIZE

	if ALL && CELLSIDFNAME == "" {
		cellIDis0 = true
	}

	for i := 0; i < limit; i++ {
		split = strings.Split(bufferarray[i], "\t")
		chro = strings.TrimPrefix(split[0], "chr")

		if cellIDis0 {
			cellID = 0
		} else {
			if cellID, isInside = celldict[split[3]];!isInside {
				continue
			}
		}

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		if _, isInside = utils.CHRINTERVALDICTTHREAD[thread][chro];!isInside {
			continue
		}

		if useOrientation {
			if len(split) <= COLSEQID {
				panic(fmt.Sprintf("Error: read split %s is out of range for orientation ID (col %d)", split, readOrientation))
			}

			readOrientation = split[COLSEQID]
		}

		if doShiftting {
			switch {
			case useOrientation && readOrientation == "-":
				start, end = start - SHIFTREAD, end - SHIFTREAD
			default:
				start, end = start + SHIFTREAD, end + SHIFTREAD
			}
		}

		intervals = utils.CHRINTERVALDICTTHREAD[thread][chro].Get(
			utils.IntInterval{Start: start, End: end})

		for _, inter = range intervals {
			if useOrientation {
				if refOrientation, isInside = ORIENTATIONDICT[uint(inter.ID())];!isInside {
					panic(fmt.Sprintf("#### ref PEAK %s (id: %d) does not have orientation!",
						utils.INTERVALMAPPING[inter.ID()], inter.ID()))
				}

				if refOrientation != readOrientation {
					continue
				}
			}

			itrg = inter.Range()

			MUTEX.Lock()

			for j = start; j < end;j++ {
				index = j - (itrg.Start + FLANKSIZE)
				switch {
				case index >= 0 && index < indexLimit:
					BASECOVERAGE[cellID][index]++
				case index < 0 && -FLANKSIZE < index:
					FLANKCOVERAGE[cellID]++
				case index > indexLimit && index < 2 * TSSREGION - FLANKSIZE:
					FLANKCOVERAGE[cellID]++
				}
			}

			if CREATEMAT {
				intervalID = inter.ID()
				initCoverageMats(intervalID)

				basecoverage = BASECOVERAGEMAT[intervalID]
				flankcoverage = FLANKCOVERAGEMAT[intervalID]

				for j = start; j < end;j++ {
					index = j - (itrg.Start + FLANKSIZE)

					switch {
					case index >= 0 && index < indexLimit:
						basecoverage[cellID][index + FLANKSIZE]++
					case index < 0 && -FLANKSIZE < index:
						flankcoverage[cellID]++
						basecoverage[cellID][index + FLANKSIZE]++
					case index > indexLimit && index < 2 * TSSREGION - FLANKSIZE:
						flankcoverage[cellID]++
						basecoverage[cellID][index + FLANKSIZE]++
					}
				}

			}

			MUTEX.Unlock()
		}
	}

	freeThreads <- thread
}


func writeCellTSSScore() {
	var buffer bytes.Buffer
	var err error
	var writer io.WriteCloser

	if STDOUT {
		writer = os.Stdout
	} else {
		writer = utils.ReturnWriter(FILENAMEOUT)
	}

	defer utils.CloseFile(writer)

	for cellID, tss := range CELLTSS {
		buffer.WriteString(cellID)
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.FormatFloat(tss, 'f', 10, 64))
		buffer.WriteRune('\n')

		if ALL {
			fmt.Printf("\n#### Global (group %s) TSS enrichment score:  %f ####\n",
				cellID, tss)
		}
	}

	_, err = writer.Write(buffer.Bytes())
	utils.Check(err)
	buffer.Reset()

	fmt.Printf("File written: %s\n", FILENAMEOUT)
}

func writeTSSMatrix() {
	waitingRef := sync.WaitGroup{}

	guardRef := make(chan bool, THREADNB)

	for i := 0;i < THREADNB;i++ {
		guardRef <- true
	}

	for clusterID := range CELLDICT {
		<- guardRef
		waitingRef.Add(1)
		go writeTSSMatrixOneGroup(clusterID, &waitingRef, guardRef)
	}

	waitingRef.Wait()
}

func writeTSSMatrixOneGroup(clusterName string, waitingRef * sync.WaitGroup, guardRef chan bool) {
	defer waitingRef.Done()

	waiting := sync.WaitGroup{}
	var intervalID uintptr

	clusterID := CELLDICT[clusterName]
	ext := path.Ext(FILENAMEOUT)
	matrixOutFile := fmt.Sprintf("%s.group_%s.deepTools.mat.gz",
		FILENAMEOUT[:len(FILENAMEOUT) - len(ext)], clusterName)

	writer := utils.ReturnWriter(matrixOutFile)
	defer utils.CloseFile(writer)

	sizeVector := (2 * TSSREGION) / int(MATRIXBINSIZE)

	header := fmt.Sprintf(`@{"verbose":false,"scale":1,"skip zeros":false,"nan after end":false,"sort using":"mean","unscaled 5 prime":[0],"body":[0],"sample_labels":["%s"],"downstream":[%d],"unscaled 3 prime":[0],"group_labels":["genes"],"bin size":[%d],"upstream":[%d],"group_boundaries":[0,%d],"sample_boundaries":[0,%d],"max threshold":null,"ref point":["center"],"min threshold":null,"sort regions":"keep","proc number":1,"bin avg type":"mean","missing data as zero":false}`, clusterName, TSSREGION, MATRIXBINSIZE, TSSREGION, len(BASECOVERAGEMAT), sizeVector)

	writer.Write([]byte(header))
	writer.Write([]byte("\n"))

	guard := make(chan bool, THREADNB)

	for i := 0;i < THREADNB;i++ {
		guard <- true
	}

	for intervalID = range BASECOVERAGEMAT {
		<- guard
		waiting.Add(1)
		// go writeOneVectorForTSSMatrix(intervalID, clusterID, &writer, &waiting, guard)
		writeOneVectorForTSSMatrix(intervalID, clusterID, &writer, &waiting, guard)
	}

	waiting.Wait()
	fmt.Printf("Matrix file written: %s\n", matrixOutFile)
	guardRef <- true
}

func writeOneVectorForTSSMatrix(intervalID uintptr,
	clusterID int,
	writer * io.WriteCloser,
	waiting * sync.WaitGroup,
	guard chan bool) {
	defer waiting.Done()

	buffer := bytes.Buffer{}
	peakstrsplit := strings.Split(utils.INTERVALMAPPING[intervalID], "\t")

	buffer.WriteString("chr")
	buffer.WriteString(strings.Join(peakstrsplit[:3], "\t"))
	buffer.WriteString("\tchr")
	buffer.WriteString(peakstrsplit[0])
	buffer.WriteRune(':')
	buffer.WriteString(peakstrsplit[1])
	buffer.WriteRune('-')
	buffer.WriteString(peakstrsplit[2])
	buffer.WriteRune('\t')
	buffer.WriteRune('.')
	buffer.WriteRune('\t')
	buffer.WriteRune('.')
	buffer.WriteRune('\t')


	sizeVector := (2 * TSSREGION) / int(MATRIXBINSIZE)

	vector := make([]string, sizeVector)

	var count uint
	var value float64
	index := 0

	var currentValue int

	basecoverage := BASECOVERAGEMAT[intervalID][clusterID]
	flankNormFactor := float64(2 * FLANKSIZE)
	flankNorm := math.Max(1, float64(FLANKCOVERAGEMAT[intervalID][clusterID]) / flankNormFactor)

	for i := 0; i < 2 * TSSREGION; i++ {
		currentValue += basecoverage[i]
		count++

		if count >= MATRIXBINSIZE {
			value = float64(currentValue) / (flankNorm * float64(count))
			vector[index] = strconv.FormatFloat(value, 'f', 2, 64)
			index++
			count = 0
			currentValue = 0
		}
	}

	if count > 0 {
		value = float64(currentValue) / (flankNorm * float64(count))
		vector[index] = strconv.FormatFloat(value, 'f', 2, 64)
	}

	buffer.WriteString(strings.Join(vector, "\t"))
	buffer.WriteRune('\n')

	MUTEX.Lock()
	_, err := (*writer).Write(buffer.Bytes())
	utils.Check(err)
	MUTEX.Unlock()

	guard <- true
}


func normaliseCountMatrix() {

	tStart := time.Now()
	chunk := len(CELLDICT) / THREADNB

	var count int
	var waiting sync.WaitGroup

	cellList := []string{}

	for cell := range CELLDICT {
		cellList = append(cellList, cell)

		count++

		if count > chunk {
			waiting.Add(1)
			go normaliseCountMatrixOneThread(cellList, &waiting)
			count = 0
			cellList = []string{}
		}
	}

	if len(cellList) > 0 {
		waiting.Add(1)
		go normaliseCountMatrixOneThread(cellList, &waiting)
	}

	waiting.Wait()

	tDiff := time.Since(tStart)
	fmt.Printf("Normalising done in time: %f s \n", tDiff.Seconds())
}


func normaliseCountMatrixOneThread(cellIDList []string, waiting *sync.WaitGroup) {
	var flankNorm, tss, maxTss float64
	var cellName string
	var i, j, cellID int

	defer waiting.Done()

	tssDict := make(map[string]float64)
	halfflank := int(SMOOTHINGWINDOW / 2)
	tssPos := TSSREGION - FLANKSIZE
	flankNormFactor := float64(2 * FLANKSIZE)

	if ALL {
		fmt.Printf("* Flank coverage: %d\n* flank size: %d\n* flank Norm factor: %f\n",
			FLANKCOVERAGE[cellID],
			int(flankNormFactor),
			float64(FLANKCOVERAGE[cellID]) / flankNormFactor)

	}

	for _, cellName = range cellIDList {
		cellID = CELLDICT[cellName]
		flankNorm = math.Max(1.0,  float64(FLANKCOVERAGE[cellID]) / flankNormFactor)
		maxTss = 0

		for i = tssPos - TSSFLANKSEARCH; i < tssPos + TSSFLANKSEARCH; i++ {
			tss = 0

			for j = i - halfflank; j < i + halfflank;j++ {
				tss += float64(BASECOVERAGE[cellID][j])
			}

			tss = tss / (float64(SMOOTHINGWINDOW))

			if tss > maxTss {
				maxTss = tss
			}
		}

		tssDict[cellName] = maxTss / (flankNorm)
	}

	if ALL {
		fmt.Printf("* Max TSS smoothed coverage: %f\n", maxTss)
	}

	MUTEX.Lock()

	for cellName, maxTss = range tssDict {
		CELLTSS[cellName] = maxTss
	}

	MUTEX.Unlock()
}
