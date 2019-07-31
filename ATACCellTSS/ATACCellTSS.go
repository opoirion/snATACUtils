/* Create excecutable to compute individual TSS for each cell */

package main


import(
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"flag"
	"os"
	"fmt"
	"log"
	"sync"
	"strings"
	"strconv"
	"github.com/biogo/store/interval"
	"time"
	"bytes"
	"path"
)


/*BASECOVERAGE map[cell][]count */
var BASECOVERAGE map[string][]int

/*FLANKCOVERAGE map[cell]count */
var FLANKCOVERAGE map[string]int

/*CELLTSS map[cell]tss */
var CELLTSS map[string]float64

/*CELLDICT map[cell]bool */
var CELLDICT map[string]bool

/*BEDFILENAME bed file name (input) */
var BEDFILENAME utils.Filename

/*PEAKFILE bed file name (input) */
var PEAKFILE utils.Filename

/*TSSFILE bed file name (input) */
var TSSFILE utils.Filename

/*CELLSIDFNAME file name file with ordered cell IDs (one ID per line) */
var CELLSIDFNAME utils.Filename

/*FILENAMEOUT  output file name output */
var FILENAMEOUT string

/*FLANKSIZE int both at the end and begining of TSS region*/
var FLANKSIZE int

/*TSSREGION int both at the end and begining of TSS region*/
var TSSREGION int

/*SMOOTHINGWINDOW int smoothing window size*/
var SMOOTHINGWINDOW int

/*TSSFLANKSEARCH tss flank size*/
var TSSFLANKSEARCH int

/*THREADNB number of threads for reading the bam file */
var THREADNB int

/*BUFFERSIZE buffer size for multithreading */
const BUFFERSIZE = 100000

/*BUFFERARRAY buffer array storing lines of bed files*/
var BUFFERARRAY[][BUFFERSIZE]string

/*MUTEX global mutex*/
var MUTEX sync.Mutex


func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `USAGE: ATACCellTSS -bed <filename> -ygi/tss <filename> -xgi <filename>
                    (optionnal -out <string> -flank <int> -smoothing <int> -boundary <int> -flank_size).

`)
		 flag.PrintDefaults()
	}

	flag.Var(&BEDFILENAME, "bed", "name of the bed file")
	flag.Var(&PEAKFILE, "ygi", "name of the bed file containing the TSS regions (assuming TSS is at the center)")
	flag.Var(&TSSFILE, "tss", "name of the bed file containing the TSS position (alternative to ygi)")
	flag.Var(&CELLSIDFNAME, "xgi", "name of the file containing the ordered list of cell IDs (one ID per line)")
	flag.StringVar(&FILENAMEOUT, "out", "", "name of the output file")
	flag.IntVar(&FLANKSIZE, "flank", 100, "flank size at the end and begining of the TSS regions")
	flag.IntVar(&TSSREGION, "boundary", 2000, "TSS boundary size at the end and begining of the TSS (used only when -tss is provided)")
	flag.IntVar(&SMOOTHINGWINDOW, "smoothing", 50, "Smoothing window size")
	flag.IntVar(&TSSFLANKSEARCH, "tss_flank", 50, "search hightest TSS values to define TSS score using this flank size arround TSS")
	flag.IntVar(&THREADNB, "threads", 1, "threads concurrency")
	flag.Parse()

	switch {
	case BEDFILENAME == "":
		log.Fatal("-bed must be provided!")
	case CELLSIDFNAME == "":
		log.Fatal("-xgi must be provided!")
	case PEAKFILE == "" && TSSFILE == "":
		log.Fatal("either -ygi or -tss must be provided!")

	case FILENAMEOUT == "":
		ext := path.Ext(CELLSIDFNAME.String())
		FILENAMEOUT = fmt.Sprintf("%s.tss_per_cell",
			CELLSIDFNAME[:len(CELLSIDFNAME)-len(ext)])

	}

	MUTEX = sync.Mutex{}

	if PEAKFILE != "" {
		utils.LoadPeaks(PEAKFILE)
	} else {
		loadTSS(TSSFILE)
	}

	initCellDicts()

	utils.CreatePeakIntervalTree()
	utils.InitIntervalDictsThreading(THREADNB)

	scanBedFile()
	normaliseCountMatrix()
	writeCellTSSScore()
}

func loadTSS(tssFile utils.Filename) {

	scanner, file := tssFile.ReturnReader(0)
	defer utils.CloseFile(file)

	var count uint
	var start int
	var buffer bytes.Buffer
	var split []string
	var err error

	utils.PEAKIDDICT = make(map[string]uint)
	count = 0

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), "\t")
		buffer.WriteString(split[0])
		buffer.WriteRune('\t')

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		buffer.WriteString(strconv.Itoa(start - TSSREGION))
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.Itoa(start + TSSREGION))
		buffer.WriteRune('\t')
		buffer.WriteString(split[3])
		buffer.WriteRune('\n')

		utils.PEAKIDDICT[buffer.String()] = count
		buffer.Reset()
		count++
	}
}

func initCellDicts() {
	var length int

	CELLDICT = utils.LoadCellIDDict(CELLSIDFNAME.String())
	FLANKCOVERAGE = make(map[string]int)
	BASECOVERAGE = make(map[string][]int)
	CELLTSS = make(map[string]float64)

	TSSREGION = assertPeakRegionHaveSameLengths()

	// if PEAKFILE != "" {
	// 	TSSREGION = assertPeakRegionHaveSameLengths()
	// }

	length = 2 * TSSREGION - 2 * FLANKSIZE

	if 2 * FLANKSIZE >= TSSREGION {
		log.Fatal(fmt.Sprintf("Error 2 * %d (flank size) greater than %d (TSS region) \n",
			FLANKSIZE, TSSREGION))
	}

	if TSSREGION <= TSSFLANKSEARCH {
		log.Fatal(fmt.Sprintf("Error TSS region (-boundary):%d must be greater than tss flank search window (-tss_flank):%d \n",
			TSSREGION, TSSFLANKSEARCH))
	}

	fmt.Printf("TSS region length: %d actual screening window size for TSS score: %d\n",
		TSSREGION, length)

	for cell := range CELLDICT {
		BASECOVERAGE[cell] = make([]int, length)
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
Error for TSS region: %s, expected to have a size of %d, found %d instead \n`,
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

func processOneBuffer(
	bufferarray *[BUFFERSIZE]string,
	thread, limit int,
	freeThreads chan int, waiting *sync.WaitGroup) {
	var split []string
	var chro, cellID string
	var start, end, j, index int
	var err error
	var isInside bool
	var intervals []interval.IntInterface
	var inter interval.IntInterface
	var itrg interval.IntRange

	defer waiting.Done()

	indexLimit := 2 * TSSREGION - 2 * FLANKSIZE

	for i := 0; i < limit; i++ {
		split = strings.Split(bufferarray[i], "\t")
		chro = split[0]
		cellID = split[3]

		if _, isInside = CELLDICT[cellID];!isInside {
			continue
		}

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		intervals = utils.CHRINTERVALDICTTHREAD[thread][chro].Get(
			utils.IntInterval{Start: start, End: end})

		for _, inter = range intervals {
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

			MUTEX.Unlock()
		}
	}

	freeThreads <- thread
}


func writeCellTSSScore() {
	var buffer bytes.Buffer
	var err error

	writer := utils.ReturnWriter(FILENAMEOUT)
	defer utils.CloseFile(writer)

	for cellID, tss := range CELLTSS {
		buffer.WriteString(cellID)
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.FormatFloat(tss, 'f', 10, 64))
		buffer.WriteRune('\n')
	}

	_, err = writer.Write(buffer.Bytes())
	utils.Check(err)
	buffer.Reset()

	fmt.Printf("File written: %s\n", FILENAMEOUT)
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

	waiting.Add(1)
	go normaliseCountMatrixOneThread(cellList, &waiting)
	waiting.Wait()

	tDiff := time.Since(tStart)
	fmt.Printf("Normalising done in time: %f s \n", tDiff.Seconds())
}


func normaliseCountMatrixOneThread(cellIDList []string, waiting *sync.WaitGroup) {
	var flankNorm, tss, maxTss float64
	var cellID string
	var i, j int

	defer waiting.Done()

	tssDict := make(map[string]float64)
	halfflank := SMOOTHINGWINDOW / 2
	tssPos := TSSREGION - FLANKSIZE

	for _, cellID = range cellIDList {
		flankNorm = float64(FLANKCOVERAGE[cellID] + 1) / 200
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

		tssDict[cellID] = maxTss / (flankNorm)
	}

	MUTEX.Lock()

	for cellID, maxTss = range tssDict {
		CELLTSS[cellID] = maxTss
	}

	MUTEX.Unlock()
}
