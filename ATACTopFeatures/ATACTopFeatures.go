package main


import (
	"flag"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"fmt"
	"log"
	"strings"
	"strconv"
	"github.com/biogo/store/interval"
	"sync"
	"time"
	stats "github.com/glycerine/golang-fisher-exact"
	"sort"
	"bytes"
)


type peak [3]string

type peakFeature struct {
	id uintptr
	cluster string
	pvalue float64
	qvalue float64
	n11, n12, n21, n22 int
}

type boolsparsematfeature struct {
	id uintptr
	cell string
}

/*BEDFILENAME bed file name (input) */
var BEDFILENAME utils.Filename

/*PEAKFILE peak file name (input) */
var PEAKFILE utils.Filename

/*CLUSTERFILE cluster file name (input) */
var CLUSTERFILE utils.Filename

/*FILENAMEOUT  output file name output */
var FILENAMEOUT string

/*CHI2ANALYSIS do chi2 analysis*/
var CHI2ANALYSIS bool

/*CONTINGENCYTABLE do chi2 analysis*/
var CONTINGENCYTABLE bool

/*THREADNB number of threads for reading the bam file */
var THREADNB int

/*THREADSCHANNEL  thread ID->channel*/
var THREADSCHANNEL chan int

/*CELLCLUSTERID  cellID <string> -> cluster ID <string>*/
var CELLCLUSTERID map[string]string

/*CLUSTERSUM  nb of cells in each cluster <string> -> sum*/
var CLUSTERSUM map[string]int

/*TOTALNBCELLS  total nb of cells*/
var TOTALNBCELLS int

/*BUFFERSIZE buffer size for multithreading */
const BUFFERSIZE = 1000000

/*SMALLBUFFERSIZE buffer size for multithreading */
const SMALLBUFFERSIZE = 500000

/*BUFFERARRAY map slice: [THREAD][BUFFERSIZE]line*/
var BUFFERARRAY [][BUFFERSIZE]string

/*BUFFERRESARRAY map slice: [THREAD][BUFFERSIZE]line*/
var BUFFERRESARRAY [][BUFFERSIZE]boolsparsematfeature

/*BOOLSPARSEMATRIX peaks x cells sparse matrix  */
var BOOLSPARSEMATRIX map[uintptr]map[string]bool

/*CHI2SCORE map[cluster ID][peak ID]count*/
var CHI2SCORE map[string][]peakFeature

/*PEAKMAPPING ...*/
var PEAKMAPPING map[uintptr]peak

/*ALPHA decision threshold*/
var ALPHA float64

/*MUTEX global mutex*/
var MUTEX sync.Mutex

/*WRITEALL Write all features even the not significant one*/
var WRITEALL bool


func main() {
	flag.Var(&BEDFILENAME, "bed", "name of the bed file")
	flag.Var(&PEAKFILE, "peak", "File containing peaks")
	flag.Var(&CLUSTERFILE, "cluster", "File containing cluster")
	flag.StringVar(&FILENAMEOUT, "out", "", "name the output file(s)")
	flag.IntVar(&THREADNB, "threads", 1, "threads concurrency")
	flag.Float64Var(&ALPHA, "alpha", 0.05, "Decision threshold")
	flag.BoolVar(&WRITEALL, "write_all", false, "Write all features including the not significant ones")
	flag.BoolVar(&CHI2ANALYSIS, "chi2", false, `perform chi2 analysis with multiple test correction
                USAGE: ATACTopFeatures -chi2 -bed <fname> -peak <fname> -cluster <fname> (optionnal -out <string> -threads <int> -alpha <float>)`)

	flag.BoolVar(&CONTINGENCYTABLE, "contingency_table", false, `Write contingency table for each feature and each cluster
                USAGE: ATACTopFeatures -chi2 -bed <fname> -peak <fname> -cluster <fname> (optionnal -out <string> -threads <int>)`)

	flag.Parse()

	switch {
	case BEDFILENAME == "":
		log.Fatal("-bed must be provided!\n")

	case PEAKFILE == "":
		log.Fatal("-peak must be provided!\n")

	case CLUSTERFILE == "":
		log.Fatal("-cluster must be provided!\n")
	}

	switch {

	case CONTINGENCYTABLE:
		createContingencyTable()
	case CHI2ANALYSIS:
		launchChi2Analysis()
	}
}


func initBoolSparseMatrix() {
	var peakKey uintptr

	BOOLSPARSEMATRIX = make(map[uintptr]map[string]bool)

	for peakKey = range PEAKMAPPING {
		BOOLSPARSEMATRIX[peakKey] = make(map[string]bool)
	}
}

func createContingencyTable() {
	BUFFERARRAY = make([][BUFFERSIZE]string, THREADNB)
	BUFFERRESARRAY = make([][BUFFERSIZE]boolsparsematfeature, THREADNB)

	if FILENAMEOUT == "" {
		FILENAMEOUT = fmt.Sprintf("%s.contingency_table.tsv\n", BEDFILENAME)
	}

	loadCellClusterIDAndInitMaps()
	utils.LoadPeaks(PEAKFILE.String())
	utils.CreatePeakIntervalTree()
	utils.InitIntervalDictsThreading(THREADNB)
	createPeakMappingDict()
	initBoolSparseMatrix()
	scanBedFile()
	computeChi2Score()
	go clearSparseMat()
	writeContingencyTable()

}


func launchChi2Analysis() {
	BUFFERARRAY = make([][BUFFERSIZE]string, THREADNB)
	BUFFERRESARRAY = make([][BUFFERSIZE]boolsparsematfeature, THREADNB)

	if FILENAMEOUT == "" {
		FILENAMEOUT = fmt.Sprintf("%s.chi2.tsv\n", BEDFILENAME)
	}

	loadCellClusterIDAndInitMaps()
	utils.LoadPeaks(PEAKFILE.String())
	utils.CreatePeakIntervalTree()
	utils.InitIntervalDictsThreading(THREADNB)
	createPeakMappingDict()
	initBoolSparseMatrix()
	scanBedFile()
	computeChi2Score()
	go clearSparseMat()
	performMultipleTestCorrection()
	writeOutput()

}


func createPeakMappingDict() {
	var peakstr string
	var peakid uint
	var peaktuple []string

	PEAKMAPPING = make(map[uintptr]peak)

	for peakstr, peakid = range utils.PEAKIDDICT {

		peaktuple = strings.Split(peakstr, "\t")
		PEAKMAPPING[uintptr(peakid)] = [3]string{
			peaktuple[0], peaktuple[1], peaktuple[2]}
	}

}


func loadCellClusterIDAndInitMaps() {
	var line, cellID, cluster string
	var split []string

	scanner, file := CLUSTERFILE.ReturnReader(0)
	defer utils.CloseFile(file)

	CELLCLUSTERID = make(map[string]string)
	CLUSTERSUM = make(map[string]int)
	CHI2SCORE = make(map[string][]peakFeature)

	for scanner.Scan() {
		line = scanner.Text()

		if line[0] == '#' {
			continue
		}

		split = strings.Split(line, "\t")
		cellID = split[0]
		cluster = split[1]

		if len(split) < 2 {
			panic(fmt.Sprintf("line: %s cannot be splitted with <tab>\n", line))
		}

		CELLCLUSTERID[cellID] = cluster
		CLUSTERSUM[cluster]++

		TOTALNBCELLS++
	}
}

func scanBedFile() {
	fmt.Printf("Scanning bed file: %s\n", BEDFILENAME)
	THREADSCHANNEL = make(chan int, THREADNB)
	tStart := time.Now()

	scanner, file := BEDFILENAME.ReturnReader(0)
	defer utils.CloseFile(file)

	count := 0

	for i:=0;i<THREADNB;i++ {
		THREADSCHANNEL <- i
	}

	threadnb := <- THREADSCHANNEL
	waiting := &sync.WaitGroup{}

	for scanner.Scan() {
		BUFFERARRAY[threadnb][count] = scanner.Text()
		count++

		if count >= BUFFERSIZE {
			waiting.Add(1)
			go processBufferArray(&BUFFERARRAY[threadnb], count, threadnb, waiting)
			count = 0
			threadnb = <- THREADSCHANNEL
		}
	}

	waiting.Add(1)
	go processBufferArray(&BUFFERARRAY[threadnb], count, threadnb, waiting)
	waiting.Wait()

	tDiff := time.Since(tStart)
	fmt.Printf("Scanning bed file done in time: %f s \n", tDiff.Seconds())
}

func processBufferArray(lineArray * [BUFFERSIZE]string, nbLines, threadnb int, waiting * sync.WaitGroup) {
	var intervals []interval.IntInterface
	var interval interval.IntInterface
	var start, end, count int
	var split []string
	var err error
	var isInside bool

	defer waiting.Done()

	for i := 0; i < nbLines; i++ {
		split = strings.Split(lineArray[i], "\t")

		if _, isInside = utils.CHRINTERVALDICT[split[0]];!isInside {
			continue
		}

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		intervals = utils.CHRINTERVALDICTTHREAD[threadnb][split[0]].Get(
			utils.IntInterval{Start: start, End: end})

		for _, interval = range intervals {
			BUFFERRESARRAY[threadnb][count].cell = split[3]
			BUFFERRESARRAY[threadnb][count].id = interval.ID()

			count++

			if count > BUFFERSIZE {
				MUTEX.Lock()
				for i :=0; i<count;i++ {
					BOOLSPARSEMATRIX[BUFFERRESARRAY[threadnb][i].id][BUFFERRESARRAY[threadnb][i].cell] = true
				}
				count = 0

				MUTEX.Unlock()

			}
		}
	}

	MUTEX.Lock()

	for i :=0; i<count;i++ {
		BOOLSPARSEMATRIX[BUFFERRESARRAY[threadnb][i].id][BUFFERRESARRAY[threadnb][i].cell] = true
	}

	MUTEX.Unlock()

	THREADSCHANNEL <- threadnb
}


func computeChi2Score() {
	fmt.Printf("Computing chi2 exact test: %s\n", BEDFILENAME)
	THREADSCHANNEL = make(chan int, THREADNB)
	var peakID uintptr

	tStart := time.Now()
	count := 0

	var chi2Array [SMALLBUFFERSIZE]uintptr

	for i:=0;i<THREADNB;i++ {
		THREADSCHANNEL <- i
	}

	threadnb := <- THREADSCHANNEL
	waiting := &sync.WaitGroup{}

	for peakID = range BOOLSPARSEMATRIX {
		chi2Array[count] = peakID

		count++

		if count > SMALLBUFFERSIZE {
			waiting.Add(1)
			go chi2ScoreOneThread(&chi2Array, waiting, count, threadnb)
			threadnb = <- THREADSCHANNEL
			count = 0
		}
	}

	waiting.Add(1)
	go chi2ScoreOneThread(&chi2Array, waiting, count, threadnb)

	waiting.Wait()

	tDiff := time.Since(tStart)
	fmt.Printf("Computing chi2 score done in time: %f s \n", tDiff.Seconds())

}


func chi2ScoreOneThread(chi2Array * [SMALLBUFFERSIZE]uintptr, waiting * sync.WaitGroup, max, threadnb int) {
	var peakID uintptr
	var cellID, cluster string
	var isInside bool
	var pvalue float64
	var value, peakTotal int
	var clusterDict map[string]int

	var results [SMALLBUFFERSIZE]peakFeature

	count := 0

	defer waiting.Done()

	for _, peakID = range (*chi2Array)[:max] {
		clusterDict = make(map[string]int)

		peakTotal = 0

		cellIDloop:
		for cellID = range BOOLSPARSEMATRIX[peakID] {
			if cluster, isInside = CELLCLUSTERID[cellID];!isInside {
				continue cellIDloop
			}

			clusterDict[cluster]++
			peakTotal++
		}

		for cluster, value = range clusterDict {
			results[count].id = peakID
			results[count].cluster = cluster

			if CONTINGENCYTABLE {
				results[count].n11 = value
				results[count].n12 = CLUSTERSUM[cluster]
				results[count].n21 = peakTotal - value
				results[count].n22 = TOTALNBCELLS - CLUSTERSUM[cluster]
			} else {

			_, pvalue = chiSquareTest(
				value, CLUSTERSUM[cluster],
				peakTotal - value, TOTALNBCELLS - CLUSTERSUM[cluster], true)

				results[count].pvalue = pvalue
			}

			count++

			if count >= BUFFERSIZE {
				MUTEX.Lock()
				for i := 0; i < count;i++ {
					CHI2SCORE[results[i].cluster] = append(CHI2SCORE[results[i].cluster], results[i])
				}

				MUTEX.Unlock()
				count = 0
			}
		}

	}

	MUTEX.Lock()
	for i := 0; i < count;i++ {
		CHI2SCORE[results[i].cluster] = append(CHI2SCORE[results[i].cluster], results[i])
	}

	MUTEX.Unlock()

	THREADSCHANNEL <- threadnb

}

func chiSquareTest(n11, n12, n21, n22 int, correction bool) (score, pvalue float64) {
	defer recoverChiScore(n11, n12, n21, n22)

	score, pvalue =  stats.ChiSquareTest(n11, n12, n21, n22, correction)

	return score, pvalue

}

func recoverChiScore(n11, n12, n21, n22 int) {

	if r := recover(); r != nil {
		fmt.Printf("ERROR: %s TABLE: %d %d %d %d\n", r, n11, n12, n21, n22)
	}
}

func performMultipleTestCorrection() {
	sortChi2Score()
	performCorrectionAfterSorting()
}

func sortChi2Score() {
	tStart := time.Now()
	THREADSCHANNEL = make(chan int, THREADNB)

	for i:=0;i<THREADNB;i++ {
		THREADSCHANNEL <- i
	}

	threadnb := <- THREADSCHANNEL

	for cluster := range CHI2SCORE {
		go func(cluster string, threadnb int) {
			sort.Slice(CHI2SCORE[cluster], func(i int, j int) bool {
				return CHI2SCORE[cluster][i].pvalue <  CHI2SCORE[cluster][j].pvalue})
			THREADSCHANNEL <- threadnb
		}(cluster, threadnb)

	threadnb = <- THREADSCHANNEL

	}

	tDiff := time.Since(tStart)
	fmt.Printf("Sorting pvalue done in time: %f s \n", tDiff.Seconds())
}

func performCorrectionAfterSorting() {

	tStart := time.Now()
	THREADSCHANNEL = make(chan int, THREADNB)

	for i:=0;i<THREADNB;i++ {
		THREADSCHANNEL <- i
	}

	threadnb := <- THREADSCHANNEL

	for cluster := range CHI2SCORE {
		go func(threadnb int, cluster string) {
			var rank int

			cllength := float64(len(CHI2SCORE[cluster]))

			for rank = range CHI2SCORE[cluster] {
				CHI2SCORE[cluster][rank].qvalue = float64(rank + 1) * ALPHA / cllength
			}

			THREADSCHANNEL <- threadnb

		}(threadnb, cluster)

		threadnb = <- THREADSCHANNEL
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Correction done in time: %f s \n", tDiff.Seconds())
}

func writeOutput() {
	var buffer bytes.Buffer
	var peakl peak
	var err error
	var isSignificant bool

	writer := utils.ReturnWriter(FILENAMEOUT)
	defer utils.CloseFile(writer)
	tStart := time.Now()

	buffer.WriteString("#chr\tstart\tstop\tcluster\tpvalue\tqvalue\tsignificant\n")

	clusters := []string{}

	for cluster := range CHI2SCORE {
		clusters = append(clusters, cluster)
	}

	sort.Strings(clusters)

	for _, cluster := range clusters {
		featureLoop:
		for _, peaki := range CHI2SCORE[cluster] {
			isSignificant = peaki.pvalue < peaki.qvalue

			if !WRITEALL && !isSignificant {
				continue featureLoop
			}

			peakl = PEAKMAPPING[peaki.id]
			buffer.WriteString(peakl[0])
			buffer.WriteRune('\t')
			buffer.WriteString(peakl[1])
			buffer.WriteRune('\t')
			buffer.WriteString(peakl[2])
			buffer.WriteRune('\t')
			buffer.WriteString(peaki.cluster)
			buffer.WriteRune('\t')
			buffer.WriteString(fmt.Sprintf("%e", peaki.pvalue))
			buffer.WriteRune('\t')
			buffer.WriteString(fmt.Sprintf("%e", peaki.qvalue))
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.FormatBool(isSignificant))
			buffer.WriteRune('\n')
		}

		_, err = writer.Write(buffer.Bytes())
		utils.Check(err)
		buffer.Reset()
	}

	tDiff := time.Since(tStart)
	fmt.Printf("File: %s written in: %f s \n", FILENAMEOUT, tDiff.Seconds())
}

func writeContingencyTable() {
	var buffer bytes.Buffer
	var peakl peak
	var err error

	writer := utils.ReturnWriter(FILENAMEOUT)
	defer utils.CloseFile(writer)
	tStart := time.Now()

	buffer.WriteString("#chr\tstart\tstop\tcluster\tn11\tn12\tn21\tn22\n")

	clusters := []string{}

	for cluster := range CHI2SCORE {
		clusters = append(clusters, cluster)
	}

	sort.Strings(clusters)

	for _, cluster := range clusters {

		for _, peaki := range CHI2SCORE[cluster] {

			peakl = PEAKMAPPING[peaki.id]
			buffer.WriteString(peakl[0])
			buffer.WriteRune('\t')
			buffer.WriteString(peakl[1])
			buffer.WriteRune('\t')
			buffer.WriteString(peakl[2])
			buffer.WriteRune('\t')
			buffer.WriteString(peaki.cluster)
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.Itoa(peaki.n11))
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.Itoa(peaki.n12))
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.Itoa(peaki.n21))
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.Itoa(peaki.n22))
			buffer.WriteRune('\n')
		}

		_, err = writer.Write(buffer.Bytes())
		utils.Check(err)
		buffer.Reset()
	}

	tDiff := time.Since(tStart)
	fmt.Printf("File: %s written in: %f s \n", FILENAMEOUT, tDiff.Seconds())
}

func clearSparseMat() {
	BOOLSPARSEMATRIX = make(map[uintptr]map[string]bool)
}
