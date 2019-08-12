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
	"path"
	"bufio"
	"os"
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

/*FEATUREPVALUEFILE Input file for MULTIPLETESTS option analysis */
var FEATUREPVALUEFILE utils.Filename

/*PEAKSYMBOLFILE peak symbol file */
var PEAKSYMBOLFILE utils.Filename

/*FILENAMEOUT  output file name output */
var FILENAMEOUT string

/*CHI2ANALYSIS do chi2 analysis*/
var CHI2ANALYSIS bool

/*CREATECONTINGENCY do chi2 analysis*/
var CREATECONTINGENCY bool

/*MULTIPLETESTS perform multiple tests correction using contingency tables file as input*/
var MULTIPLETESTS bool

/*THREADNB number of threads for reading the bam file */
var THREADNB int

/*SPLIT split the initial list of peaks into sets for memory efficiency */
var SPLIT int

/*THREADSCHANNEL  thread ID->channel*/
var THREADSCHANNEL chan int

/*CELLCLUSTERID  cellID <string> -> cluster ID <string>*/
var CELLCLUSTERID map[string]string

/*CLUSTERSUM  nb of cells in each cluster <string> -> sum*/
var CLUSTERSUM map[string]int

/*TOTALNBCELLS  total nb of cells*/
var TOTALNBCELLS int

/*BUFFERSIZE buffer size for multithreading */
const BUFFERSIZE = 50000

/*SMALLBUFFERSIZE buffer size for multithreading */
const SMALLBUFFERSIZE = 50000

/*BUFFERARRAY map slice: [THREAD][BUFFERSIZE]line*/
var BUFFERARRAY [][BUFFERSIZE]string

/*BUFFERRESARRAY map slice: [THREAD][BUFFERSIZE]line*/
var BUFFERRESARRAY [][BUFFERSIZE]boolsparsematfeature

/*BOOLSPARSEMATRIX peaks x cells sparse matrix  */
var BOOLSPARSEMATRIX []map[string]bool

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

/*PEAKSYMBOLDICT map[peak]symbol */
var PEAKSYMBOLDICT map[peak]string


func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `
#################### MODULE TO INFER SIGNIFICANT CLUSTER PEAKS ########################

"""full individual chi2 computation for each peak with FDR correction using Benjamini-Hochberg correction. Not recommended because using golang suboptimal chi2 implementation"""
USAGE: ATACTopFeatures -chi2 -bed <fname> -peak <fname> -cluster <fname> (optionnal -out <string> -threads <int> -alpha <float> -write_all -split <int>)

"""Create contingency table for each feature and each cluster"""
USAGE: ATACTopFeatures -create_contingency -bed <fname> -peak <fname> -cluster <fname> (optionnal -out <string> -threads <int>)

"""correct feature pvalue for multiple tests performed or each cluster"""
USAGE: ATACTopFeatures -pvalue_correction -ptable <fname> (optionnal -out <string> -threads <int> -alpha <float> -write_all)


`)
		 flag.PrintDefaults()
	}


	flag.Var(&BEDFILENAME, "bed", "name of the bed file")
	flag.Var(&PEAKFILE, "peak", "File containing peaks")
	flag.Var(&PEAKSYMBOLFILE, "symbol", `File containing symbols (such as gene name) for peak file.
     Each row should either contain one symbol per line and matches the peaks from -peak OR option2:<symbol>\t<chromosome>\t<start>\t<stop>\n`)
	flag.Var(&CLUSTERFILE, "cluster", "File containing cluster")
	flag.Var(&FEATUREPVALUEFILE, "ptable", `File containing pvalue for each interval feature
                row scheme: <chromosome>\t<start>\t<stop>\t<cluster ID>\t<pvalue>\n`)
	flag.StringVar(&FILENAMEOUT, "out", "", "name the output file(s)")
	flag.IntVar(&THREADNB, "threads", 1, "threads concurrency")
	flag.IntVar(&SPLIT, "split", 0, "Split the input set of peaks into multiple subsets (The number is defined by the -split option) processed one by one for memory efficiency.")
	flag.Float64Var(&ALPHA, "alpha", 0.05, "Decision threshold")
	flag.BoolVar(&WRITEALL, "write_all", false, "Write all features including the not significant ones")
	flag.BoolVar(&CHI2ANALYSIS, "chi2", false, `perform chi2 analysis with multiple test correction`)
	flag.BoolVar(&CREATECONTINGENCY, "create_contingency", false, `Create contingency table for each feature and each cluster`)
	flag.BoolVar(&MULTIPLETESTS, "pvalue_correction", false, `correct feature pvalue for multiple tests performed or each cluster`)

	flag.Parse()

	switch {
	case MULTIPLETESTS:
		if FEATUREPVALUEFILE == "" {
			log.Fatal("-ptable must be provided!\n")
		}

		launchMultipleTestAnalysis()
		return
	}

	switch {
	case BEDFILENAME == "":
		log.Fatal("-bed must be provided!\n")

	case PEAKFILE == "":
		log.Fatal("-peak must be provided!\n")

	case CLUSTERFILE == "":
		log.Fatal("-cluster must be provided!\n")
	}

	switch {

	case CREATECONTINGENCY:
		if SPLIT > 1 {
			createContingencyTableUsingSubsets()
		} else {
			createContingencyTable()
		}
	case CHI2ANALYSIS:
		launchChi2Analysis()
	}
}


func launchMultipleTestAnalysis() {
	if FILENAMEOUT == "" {
		ext := path.Ext(FEATUREPVALUEFILE.String())
		FILENAMEOUT = fmt.Sprintf("%s.pvalue_corrected.tsv",
			FEATUREPVALUEFILE[:len(FEATUREPVALUEFILE)-len(ext)])
	}

	loadSymbolFile()
	loadPvalueTable()
	sortPvalueScore()
	performCorrectionAfterSorting()
	writePvalueCorrectedTable()
}


func createContingencyTable() {
	BUFFERARRAY = make([][BUFFERSIZE]string, THREADNB)
	BUFFERRESARRAY = make([][BUFFERSIZE]boolsparsematfeature, THREADNB)

	if FILENAMEOUT == "" {
		ext := path.Ext(BEDFILENAME.String())
		FILENAMEOUT = fmt.Sprintf("%s.contingency_table.tsv",
			BEDFILENAME[:len(BEDFILENAME)-len(ext)])
	}

	tStart := time.Now()

	loadSymbolFile()
	loadCellClusterIDAndInitMaps()
	utils.LoadPeaks(PEAKFILE)
	utils.CreatePeakIntervalTree()
	utils.InitIntervalDictsThreading(THREADNB)
	createPeakMappingDict()
	initBoolSparseMatrix()
	scanBedFile()
	computeChi2Score()
	go clearSparseMat()
	writeContingencyTable(FILENAMEOUT, true)

	tDiff := time.Since(tStart)
	fmt.Printf("Create contingency table done in time: %f s \n", tDiff.Seconds())

}

func createContingencyTableUsingSubsets() {
	var chunk, firstPeak, lastPeak int
	var filenameout string
	var filenames []string

	BUFFERARRAY = make([][BUFFERSIZE]string, THREADNB)
	BUFFERRESARRAY = make([][BUFFERSIZE]boolsparsematfeature, THREADNB)

	if FILENAMEOUT == "" {
		ext := path.Ext(BEDFILENAME.String())
		FILENAMEOUT = fmt.Sprintf("%s.contingency_table.tsv",
			BEDFILENAME[:len(BEDFILENAME)-len(ext)])
	}

	loadSymbolFile()
	utils.LoadPeaks(PEAKFILE)

	chunk = (len(utils.PEAKIDDICT) + SPLIT) / SPLIT
	lastPeak = chunk

	tStart := time.Now()

	for i := 0; i < SPLIT; i ++ {
		filenameout = fmt.Sprintf("%s.%d", FILENAMEOUT, i)

		fmt.Printf("Processing set of peaks: (start: %d  end:%d)\n", firstPeak, lastPeak)
		utils.LoadPeaksSubset(PEAKFILE, firstPeak, lastPeak)
		loadSymbolFile()
		loadCellClusterIDAndInitMaps()

		utils.CreatePeakIntervalTree()
		utils.InitIntervalDictsThreading(THREADNB)

		createPeakMappingDict()
		initBoolSparseMatrix()
		scanBedFile()
		computeChi2Score()
		go clearSparseMat()
		writeContingencyTable(filenameout, i == 0)

		lastPeak += chunk
		firstPeak += chunk

		filenames = append(filenames, filenameout)
	}

	fmt.Printf("Combining contingency table and cleaning...\n")
	mergeAndCleanTmpContingencyTables(filenames)

	tDiff := time.Since(tStart)
	fmt.Printf("Create contingency table done in time: %f s \n", tDiff.Seconds())
}


func mergeAndCleanTmpContingencyTables(filenames []string) {
	var cmd, cmdRm bytes.Buffer

	cmd.WriteString("cat ")
	cmdRm.WriteString("rm ")

	for _, filename := range filenames {
		cmd.WriteString(filename)
		cmd.WriteRune(' ')

		cmdRm.WriteString(filename)
		cmdRm.WriteRune(' ')
	}

	cmd.WriteString(" > ")
	cmd.WriteString(FILENAMEOUT)

	utils.ExceCmd(cmd.String())
	utils.ExceCmd(cmdRm.String())
}

func launchChi2Analysis() {
	BUFFERARRAY = make([][BUFFERSIZE]string, THREADNB)
	BUFFERRESARRAY = make([][BUFFERSIZE]boolsparsematfeature, THREADNB)

	if FILENAMEOUT == "" {
		ext := path.Ext(BEDFILENAME.String())
		FILENAMEOUT = fmt.Sprintf("%s.chi2.tsv",
			BEDFILENAME[:len(BEDFILENAME)-len(ext)])
	}

	loadSymbolFile()
	loadCellClusterIDAndInitMaps()
	utils.LoadPeaks(PEAKFILE)
	utils.CreatePeakIntervalTree()
	utils.InitIntervalDictsThreading(THREADNB)
	createPeakMappingDict()
	initBoolSparseMatrix()
	scanBedFile()
	computeChi2Score()
	go clearSparseMat()
	performMultipleTestCorrection()
	writePvalueCorrectedTable()

}


func loadSymbolFile() {
	var scannerPeak *bufio.Scanner
	var filePeak *os.File
	var split, split2 []string
	var peakl peak
	var symbol string

	PEAKSYMBOLDICT = make(map[peak]string)

	if PEAKSYMBOLFILE == "" {
		return
	}

	isOption1 := true

	if PEAKFILE == "" {
		isOption1 = false
	} else {
		scannerPeak, filePeak = PEAKFILE.ReturnReader(0)
		defer utils.CloseFile(filePeak)
	}

	scanner, file := PEAKSYMBOLFILE.ReturnReader(0)
	defer utils.CloseFile(file)

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), "\t")

		if len(split) == 4 {
			isOption1 = false
		}

		if !isOption1 && len(split) != 4 {
			panic(fmt.Sprintf(
				"Error line %s from symbol file should be <symbol>\t<chromosome>\t<start>\t<stop>\n",
				split))
		}

		symbol = split[0]

		if isOption1 {
			scannerPeak.Scan()
			split2 = strings.Split(scannerPeak.Text(), "\t")
			peakl = peak{split2[0], split2[1], split2[2]}

		} else {
			peakl = peak{split[1], split2[2], split2[3]}
		}

		PEAKSYMBOLDICT[peakl] = symbol
	}
}


func loadPvalueTable() {
	var line, symbol string
	var split []string
	var peaki peakFeature
	var peakl peak
	var err error
	var isInside bool
	var count uintptr
	var pvalue float64

	tStart := time.Now()

	CHI2SCORE = make(map[string][]peakFeature)
	PEAKMAPPING = make(map[uintptr]peak)

	peakset := make(map[peak]uintptr)

	scanner, file := FEATUREPVALUEFILE.ReturnReader(0)
	defer utils.CloseFile(file)

	isSymbolFile := PEAKSYMBOLFILE != ""

	for scanner.Scan() {
		line = scanner.Text()
		if line[0] == '#' {
			continue
		}

		split = strings.Split(line, "\t")

		_, err = strconv.Atoi(split[1])
		utils.Check(err)

		_, err = strconv.Atoi(split[2])
		utils.Check(err)

		peakl = peak{split[0], split[1], split[2]}

		if _, isInside = peakset[peakl];!isInside {
			peakset[peakl] = count
			PEAKMAPPING[count] = peakl
			count++
		}

		peaki.id  = peakset[peakl]
		peaki.cluster = split[3]
		pvalue, err = strconv.ParseFloat(split[len(split) - 1], 64)
		utils.Check(err)
		peaki.pvalue = pvalue
		CHI2SCORE[peaki.cluster] = append(CHI2SCORE[peaki.cluster], peaki)

		if !isSymbolFile && len(split) == 6 {
			symbol = split[4]
			PEAKSYMBOLDICT[peakl] = symbol
		}
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Table %s loaded in: %f s \n", FEATUREPVALUEFILE, tDiff.Seconds())
}

func initBoolSparseMatrix() {
	var peakKey uintptr

	BOOLSPARSEMATRIX = make([]map[string]bool, len(PEAKMAPPING))

	for peakKey = range PEAKMAPPING {
		BOOLSPARSEMATRIX[peakKey] = make(map[string]bool)
	}
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
	var inter interval.IntInterface
	var start, end, count int
	var split []string
	var err error
	var isInside bool
	var intree *interval.IntTree

	defer waiting.Done()

	for i := 0; i < nbLines; i++ {
		split = strings.Split(lineArray[i], "\t")

		if _, isInside = CELLCLUSTERID[split[3]];!isInside {
			continue
		}

		if intree, isInside = utils.CHRINTERVALDICT[split[0]];!isInside {
			continue
		}

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		intervals = intree.Get(
			utils.IntInterval{Start: start, End: end})

		for _, inter = range intervals {
			BUFFERRESARRAY[threadnb][count].cell = split[3]
			BUFFERRESARRAY[threadnb][count].id = inter.ID()

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
	fmt.Printf("Computing chi2 exact test...\n")
	THREADSCHANNEL = make(chan int, THREADNB)
	var peakID int

	tStart := time.Now()
	count := 0

	chunk := SMALLBUFFERSIZE

	chi2Array := make([][]uintptr, THREADNB)

	for i:=0;i<THREADNB;i++ {
		THREADSCHANNEL <- i
		chi2Array[i] = make([]uintptr, chunk)
	}

	threadnb := <- THREADSCHANNEL
	var waiting sync.WaitGroup

	for peakID = range BOOLSPARSEMATRIX {
		chi2Array[threadnb][count] = uintptr(peakID)

		count++

		if count >= chunk {
			waiting.Add(1)
			go chi2ScoreOneThread(&chi2Array[threadnb], &waiting, count, threadnb)
			threadnb = <- THREADSCHANNEL
			count = 0
		}
	}

	waiting.Add(1)
	go chi2ScoreOneThread(&chi2Array[threadnb], &waiting, count, threadnb)

	waiting.Wait()

	tDiff := time.Since(tStart)
	fmt.Printf("Computing chi2 score done in time: %f s \n", tDiff.Seconds())

}


func chi2ScoreOneThread(chi2Array * []uintptr, waiting * sync.WaitGroup, max, threadnb int) {
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

			if CREATECONTINGENCY {
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

			if count >= SMALLBUFFERSIZE {
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
	sortPvalueScore()
	performCorrectionAfterSorting()
}

func sortPvalueScore() {
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

func writePvalueCorrectedTable() {
	var buffer bytes.Buffer
	var peakl peak
	var err error
	var isSignificant bool

	writeSymbol := len(PEAKSYMBOLDICT) != 0

	writer := utils.ReturnWriter(FILENAMEOUT)
	defer utils.CloseFile(writer)
	tStart := time.Now()

	buffer.WriteString("#chr\tstart\tstop\tcluster\tpvalue\tqvalue\tsignificant")

	if writeSymbol {
		buffer.WriteString("\tsymbol")
	}

	buffer.WriteRune('\n')

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

			if writeSymbol {
				buffer.WriteRune('\t')
				buffer.WriteString(PEAKSYMBOLDICT[peakl])
			}

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

func writeContingencyTable(filenameout string, header bool) {
	var buffer bytes.Buffer
	var peakl peak
	var err error

	writer := utils.ReturnWriter(filenameout)
	defer utils.CloseFile(writer)
	tStart := time.Now()

	writeSymbol := len(PEAKSYMBOLDICT) != 0

	if header {
		buffer.WriteString("#chr\tstart\tstop\tcluster")

		if writeSymbol {
			buffer.WriteString("\tsymbol")
		}

		buffer.WriteString("\tn11\tn12\tn21\tn22\n")
	}

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

			if writeSymbol {
				buffer.WriteRune('\t')
				buffer.WriteString(PEAKSYMBOLDICT[peakl])
			}

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
	BOOLSPARSEMATRIX = make([]map[string]bool, 0)
}
