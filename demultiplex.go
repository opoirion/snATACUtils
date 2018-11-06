package main

import "fmt"
import "flag"
import "os"
import "log"
import "bufio"
import "strings"
import "time"
import "errors"
import "sync"
import "io"
import "path"


import (
	utils "ATACdemultiplex/ATACdemultiplexUtils"
	pathutils "path"
)

/*PRINTVERSION ...*/
var PRINTVERSION bool
/*FASTQ_R1 ...*/
var FASTQ_R1 string
/*FASTQ_R2 ...*/
var FASTQ_R2 string
/*FASTQ_I1 ...*/
var FASTQ_I1 string
/*FASTQ_I2 ...*/
var FASTQ_I2 string
/*NB_THREADS ...*/
var NB_THREADS int
/*TAGLENGTH ...*/
var TAGLENGTH int
/*MAX_NB_READS ...*/
var MAX_NB_READS int
/*COMPRESSION_MODE ...*/
var COMPRESSION_MODE int
/*USE_BZIP_GO_LIBRARY ...*/
var USE_BZIP_GO_LIBRARY bool
/*WRITE_LOGS ...*/
var WRITE_LOGS bool
/*OUTPUT_TAG_NAME ...*/
var OUTPUT_TAG_NAME string
/*INDEX_REPLICATE_R1 ...*/
var INDEX_REPLICATE_R1 string
/*INDEX_REPLICATE_R2 ...*/
var INDEX_REPLICATE_R2 string
/*INDEX_NO_REPLICATE ...*/
var INDEX_NO_REPLICATE string
/*OUTPUT_PATH ...*/
var OUTPUT_PATH string
/*MAX_NB_MISTAKE ...*/
var MAX_NB_MISTAKE int
/*DEBUG ...*/
var DEBUG bool
/*ERRORHANDLING ...*/
var ERRORHANDLING string

/*INDEX_R1_DICT ...*/
var INDEX_R1_DICT map[string]map[string]bool
/*INDEX_R2_DICT ...*/
var INDEX_R2_DICT map[string]map[string]bool
/*INDEX_NO_DICT ...*/
var INDEX_NO_DICT map[string]map[string]bool


/*REPLNUMBER ...*/
var REPLNUMBER int

/*LOG_CHAN ...*/
var LOG_CHAN map[string]chan StatsDict
/*LOG_INDEX_READ_CHAN ...*/
var LOG_INDEX_READ_CHAN map[string]chan StatsDict
/*LOG_INDEX_CELL_CHAN ...*/
var LOG_INDEX_CELL_CHAN map[string]chan StatsDict
/*SORT_LOGS ...*/
var SORT_LOGS bool
/*SHIFT_P5 ...*/
var SHIFT_P5 int

/*LENGTHDIC length of the different indexes*/
var LENGTHDIC = map[string]int {
	"i5":0,
	"i7":0,
	"p5":0,
	"p7":0,
}

/*LOG_TYPE ...*/
var LOG_TYPE  = []string {
	"stats",
	"fail",
}

/*LOG_INDEX_TYPE ...*/
var LOG_INDEX_TYPE  = []string {
	"fail_p5",
	"fail_p7",
	"fail_i5",
	"fail_i7",
}



/*StatsDict doc */
type StatsDict struct {
	dict map[string]int
	name string
}

/*Init doc */
func (statDict *StatsDict) Init() {
	statDict.dict = make(map[string]int)
}


/* */
func initChan(logChan * map[string]chan StatsDict, logType []string) {

	(*logChan) = make(map[string] chan StatsDict)

	for _, value := range logType {
		(*logChan)[value] = make(chan StatsDict, NB_THREADS)
	}
}


/* */
func main() {
	// counter := make(map[string]int)
	flag.StringVar(&FASTQ_I1, "fastq_I1", "", "fastq index file index paired read 1")
	flag.StringVar(&FASTQ_I2, "fastq_I2", "", "fastq index file index paired read 2")
	flag.StringVar(&FASTQ_R1, "fastq_R1", "", "fastq read file index paired read 1")
	flag.StringVar(&FASTQ_R2, "fastq_R2", "", "fastq read file index paired read 2")
	flag.BoolVar(&DEBUG, "debug", false, "debug wrongly formated reads")
	flag.BoolVar(&PRINTVERSION, "version", false, "print the current version and return")
	flag.IntVar(&MAX_NB_MISTAKE, "max_nb_mistake", 2, "Maximum number of mistakes allowed to assign a reference read id (default 2)")
	flag.StringVar(&OUTPUT_TAG_NAME, "output_tag_name", "", "tag for the output file names (default None)")
	flag.BoolVar(&USE_BZIP_GO_LIBRARY, "use_bzip2_go_lib", false, "use bzip2 go library instead of native C lib (slower)")
	flag.BoolVar(&WRITE_LOGS, "write_logs", false, "write logs (might slower the execution time)")
	flag.BoolVar(&SORT_LOGS, "sort_logs", false, "sort logs (might consume a lot of RAM and provoke failure)")
	flag.IntVar(&SHIFT_P5, "shift_p5", 0, "shift p5 barcodes toward n nucleotides from the left (default 0)")

	flag.IntVar(&COMPRESSION_MODE, "compressionMode", 6, `compressionMode for native bzip2 lib
 (1 faster -> 9 smaller) <default: 6>`)
	flag.IntVar(&NB_THREADS, "nbThreads", 1, "number of threads to use")
	flag.IntVar(&TAGLENGTH, "taglength", 8,
		`<OPTIONAL> number of nucleotides to consider at the end
 and begining (default 8)`)
	flag.IntVar(&MAX_NB_READS, "max_nb_reads", 0,
		"<OPTIONAL> max number of reads to process (default 0 => None)")
	flag.StringVar(&INDEX_REPLICATE_R1, "index_replicate_r1", "",
		"<OPTIONAL> path toward indexes of R1 replicates (i.e. replicate number 1)")
	flag.StringVar(&INDEX_REPLICATE_R2, "index_replicate_r2", "",
		"<OPTIONAL> path toward indexes of R2 replicates (i.e. replicate number 2)")
	flag.StringVar(&INDEX_NO_REPLICATE, "index_no_replicate", "",
		"<OPTIONAL> path toward indexes when only 1 replicate is used")
	flag.StringVar(&OUTPUT_PATH, "output_path", "",
		"<OPTIONAL> output path being used")
	flag.StringVar(&ERRORHANDLING, "error_handling", "return",
		"error handling strategy (currently return or raise)." +
		" if return, the demultiplex returns in case of an error and continue.")
	flag.Parse()

	if PRINTVERSION {
		fmt.Printf("ATACdemultiplex version: %s\n", VERSION)
		return
	}

	fmt.Printf("#### current ATACdemultiplex version: %s\n", VERSION)

	switch {
	case INDEX_NO_REPLICATE != "":
		loadIndexes(INDEX_NO_REPLICATE, &INDEX_NO_DICT)
		REPLNUMBER = 1

		if INDEX_REPLICATE_R1 != "" || INDEX_REPLICATE_R2 != "" {
			log.Fatal("Cannot set up index_no_replicate with either index_replicate_r1 or index_replicate_r2")
		}

	default:
		loadIndexes(INDEX_REPLICATE_R1, &INDEX_R1_DICT)
		loadIndexes(INDEX_REPLICATE_R2, &INDEX_R2_DICT)
		REPLNUMBER = 2

		if  INDEX_NO_REPLICATE != "" {
			log.Fatal("Cannot set up index_no_replicate with either index_replicate_r1 or index_replicate_r2")
		}
		if INDEX_REPLICATE_R1 == "" || INDEX_REPLICATE_R2 == "" {
			log.Fatal("Both index_replicate_r1 and index_replicate_r2 should be set up!")
		}
	}

	for repl := 1; repl < REPLNUMBER + 1;repl++ {
		LOG_INDEX_TYPE = append(LOG_INDEX_TYPE, fmt.Sprintf("success_p5_repl%d", repl))
		LOG_INDEX_TYPE = append(LOG_INDEX_TYPE, fmt.Sprintf("success_p7_repl%d", repl))
		LOG_INDEX_TYPE = append(LOG_INDEX_TYPE, fmt.Sprintf("success_i5_repl%d", repl))
		LOG_INDEX_TYPE = append(LOG_INDEX_TYPE, fmt.Sprintf("success_i7_repl%d", repl))
		LOG_TYPE = append(LOG_TYPE, fmt.Sprintf("success_repl%d", repl))
	}

	initChan(&LOG_CHAN, LOG_TYPE)
	initChan(&LOG_INDEX_READ_CHAN, LOG_INDEX_TYPE)
	initChan(&LOG_INDEX_CELL_CHAN, LOG_INDEX_TYPE)

	if len(OUTPUT_TAG_NAME) > 0 {
		OUTPUT_TAG_NAME = fmt.Sprintf("_%s", OUTPUT_TAG_NAME)
	}

	if OUTPUT_PATH == "." || OUTPUT_PATH == "" {
		OUTPUT_PATH = "./"
	}

	OUTPUT_PATH = fmt.Sprintf("%s/", OUTPUT_PATH)

	if _, err := os.Stat(OUTPUT_PATH); err != nil {
		os.Mkdir(OUTPUT_PATH, os.ModePerm)
	}

	fmt.Printf("fastq file index 1 analyzed: %s\n", FASTQ_I1)
	fmt.Printf("fastq file index 2 analyzed: %s\n", FASTQ_I2)
	fmt.Printf("fastq file read file 1 analyzed: %s\n", FASTQ_R1)
	fmt.Printf("fastq file read file 2 analyzed: %s\n", FASTQ_R2)

	tStart := time.Now()

	switch {
	case NB_THREADS == 1:
		var waiting sync.WaitGroup
		waiting.Add(1)
		launchAnalysisOneFile(0, MAX_NB_READS, "", &waiting)

	case NB_THREADS < 1:
		panic(errors.New("threads should be >= 1! "))

	case NB_THREADS > 1:
		launchAnalysisMultipleFile()
	}

	writeReport()

	tDiff := time.Now().Sub(tStart)
	fmt.Printf("demultiplexing finished in %f s\n", tDiff.Seconds())
}

/* */
func writeReport() {

	if !WRITE_LOGS{
		return
	}
	tStart := time.Now()
	fmt.Printf("writing reports....\n")
	writeReportFromMultipleDict(&LOG_INDEX_CELL_CHAN, "index_cell")
	writeReportFromMultipleDict(&LOG_INDEX_READ_CHAN, "index_read")

	for key, logChan := range LOG_CHAN {
		logs := extractDictFromChan(logChan)
		filename := fmt.Sprintf("%sreport%s_%s.log",
			OUTPUT_PATH, OUTPUT_TAG_NAME, key)

		file, err := os.Create(filename)
		defer file.Close()
		Check(err)

		file.WriteString("#<key>\t<value>\n")

		switch SORT_LOGS {
		case true:
			rankedLogs := utils.RankByWordCountAndDeleteOldMap(&logs)

			for _, pair := range rankedLogs {
				file.WriteString(fmt.Sprintf("%s\t%d\n", pair.Key, pair.Value))
			}

		default:

			for k, v := range logs {
				file.WriteString(fmt.Sprintf("%s\t%d\n", k, v))
			}
		}
	}

	tDiff := time.Now().Sub(tStart)
	fmt.Printf("writing report finished in: %f s\n", tDiff.Seconds())
}


/* */
func writeReportFromMultipleDict(channel * map[string]chan StatsDict, fname string) {
	filename := fmt.Sprintf("%sreport%s_%s.log",
		OUTPUT_PATH, OUTPUT_TAG_NAME, fname)
	file, err := os.Create(filename)

	defer file.Close()
	Check(err)


	filename = fmt.Sprintf("%sreport%s_%s_fail.log",
		OUTPUT_PATH, OUTPUT_TAG_NAME, fname)
	fileFail, err := os.Create(filename)

	defer fileFail.Close()
	Check(err)

	var fp *os.File

	for dictType, logChan := range *channel {

		switch{
		case strings.Contains(dictType, "fail"):
			fp = fileFail

			if SORT_LOGS {
				fp.WriteString(fmt.Sprintf("#### %s\n", dictType))
			}
		default:
			fp = file
			fp.WriteString(fmt.Sprintf("#### %s\n", dictType))
		}

		logs := extractDictFromChan(logChan)

		switch SORT_LOGS {
		case true:
			rankedLogs := utils.RankByWordCountAndDeleteOldMap(&logs)

			for _, pair := range rankedLogs {
				fp.WriteString(fmt.Sprintf("%s\t%s\t%d\n", dictType, pair.Key, pair.Value))
			}

		default:

			for k, v := range logs {
				fp.WriteString(fmt.Sprintf("%s\t%s\t%d\n", dictType, k, v))
			}
		}
	}
}

/* */
func extractDictFromChan(channel chan StatsDict) (map[string]int) {

	dict := make(map[string]int)

	loop:
	for {
		select{
		case statsDict := <-channel:
			dictit := statsDict.dict
			for key, value := range dictit {
				dict[key] += value
			}
		default:
			break loop
		}
	}

	return dict
}

/* */
func launchAnalysisMultipleFile() {

	var nbReads int
	var chunk int

	switch {
	case MAX_NB_READS != 0:
		chunk = (MAX_NB_READS / NB_THREADS)
	default:
		nbReads = countLine(FASTQ_I1, COMPRESSION_MODE) / 4
		fmt.Printf("estimated number of reads: %d\n", nbReads)
		chunk = (nbReads / NB_THREADS) - (nbReads / NB_THREADS) % 4
	}

	fmt.Printf("chunck size: %d\n", chunk)

	startingRead := 0
	index := 0

	var waiting sync.WaitGroup
	waiting.Add(NB_THREADS)

	for i := 0; i < NB_THREADS; i++ {
		if i == NB_THREADS-1 && MAX_NB_READS == 0 {
			chunk = 0
		}

		go launchAnalysisOneFile(
			startingRead, chunk,
			fmt.Sprintf("index_%d.", index),
			&waiting)

		startingRead += chunk
		index++
	}

	waiting.Wait()

	_, filenameR1 := pathutils.Split(FASTQ_R1)
	_, filenameR2 := pathutils.Split(FASTQ_R2)

	ext := path.Ext(filenameR1)

	for repl := 1;repl < REPLNUMBER +1; repl++ {

		outputR1 := fmt.Sprintf("%s%s%s.demultiplexed.R1.repl%d.fastq%s", OUTPUT_PATH,
			strings.TrimSuffix(filenameR1, fmt.Sprintf(".fastq%s", ext)),
			OUTPUT_TAG_NAME, repl, ext)
		outputR2 := fmt.Sprintf("%s%s%s.demultiplexed.R2.repl%d.fastq%s", OUTPUT_PATH,
			strings.TrimSuffix(filenameR2, fmt.Sprintf(".fastq%s", ext)),
			OUTPUT_TAG_NAME, repl, ext)

		cmd1 := fmt.Sprintf("cat %sindex_*.demultiplexed.R1.repl%d.fastq%s > %s",
			OUTPUT_PATH, repl, ext, outputR1)
		cmd2 := fmt.Sprintf("cat %sindex_*.demultiplexed.R2.repl%d.fastq%s > %s",
			OUTPUT_PATH, repl, ext, outputR2)

		fmt.Printf("concatenating repl. %d read 1 files...\n", repl)
		fmt.Printf("%s\n", cmd1)
		utils.ExceCmd(cmd1)
		fmt.Printf("concatenating repl. %d read 2 files...\n", repl)
		fmt.Printf("%s\n", cmd2)
		utils.ExceCmd(cmd2)

	}

	cmd := fmt.Sprintf("rm %sindex_*demultiplexed*repl*%s ", OUTPUT_PATH, ext)

	fmt.Printf("removing read index files...\n")
	utils.ExceCmd(cmd)
}

/* */
func initLog(logType []string) (logs map[string]*StatsDict) {

	logs = make(map[string]*StatsDict)

	for _, value := range logType {
			logs[value] = &StatsDict{name:value}
			logs[value].Init()
	}

	return logs
}

/* */
func launchAnalysisOneFile(
	startingRead int,
	maxNbReads int,
	index string,
	waiting *sync.WaitGroup) {

	defer waiting.Done()

	logs := initLog(LOG_TYPE)

	logsIndexRead := initLog(LOG_INDEX_TYPE)
	logsIndexCell := initLog(LOG_INDEX_TYPE)

	var scannerI1 * bufio.Scanner
	var scannerI2 * bufio.Scanner
	var scannerR1 * bufio.Scanner
	var scannerR2 * bufio.Scanner

	var fileI1 * os.File
	var fileI2 * os.File
	var fileR1 * os.File
	var fileR2 * os.File

	_, filenameR1 := pathutils.Split(FASTQ_R1)
	_, filenameR2 := pathutils.Split(FASTQ_R2)

	ext := path.Ext(filenameR1)

	scannerI1, fileI1 = utils.ReturnReader(FASTQ_I1, startingRead * 4, USE_BZIP_GO_LIBRARY)
	scannerI2, fileI2 = utils.ReturnReader(FASTQ_I2, startingRead * 4, USE_BZIP_GO_LIBRARY)
	scannerR1, fileR1 = utils.ReturnReader(FASTQ_R1, startingRead * 4, USE_BZIP_GO_LIBRARY)
	scannerR2, fileR2 = utils.ReturnReader(FASTQ_R2, startingRead * 4, USE_BZIP_GO_LIBRARY)


	writerR1DictName := make(map[int]string)
	writerR2DictName := make(map[int]string)

	writerR1Dict := make(map[string]io.WriteCloser)
	writerR2Dict := make(map[string]io.WriteCloser)

	for repl := 1 ; repl < REPLNUMBER +1; repl ++ {
		writerR1DictName[repl] = fmt.Sprintf("%s%s%s%s.demultiplexed.R1.repl%d.fastq%s",
			OUTPUT_PATH, index,
			strings.TrimSuffix(filenameR1, fmt.Sprintf(".fastq%s", ext)),
			OUTPUT_TAG_NAME, repl, ext)
		writerR2DictName[repl] = fmt.Sprintf("%s%s%s%s.demultiplexed.R2.repl%d.fastq%s", OUTPUT_PATH, index,
			strings.TrimSuffix(filenameR2, fmt.Sprintf(".fastq%s", ext)),
			OUTPUT_TAG_NAME, repl, ext)
		writerR1Dict[writerR1DictName[repl]] = utils.ReturnWriter(
			writerR1DictName[repl], COMPRESSION_MODE,
			USE_BZIP_GO_LIBRARY)
		writerR2Dict[writerR2DictName[repl]] = utils.ReturnWriter(
			writerR2DictName[repl], COMPRESSION_MODE,
			USE_BZIP_GO_LIBRARY)

		defer writerR1Dict[writerR1DictName[repl]].Close()
		defer writerR2Dict[writerR2DictName[repl]].Close()
	}

	defer fileI1.Close()
	defer fileI2.Close()
	defer fileR1.Close()
	defer fileR2.Close()

	var id_I1, id_I2, id_R1, id_R2 string
	var read_I1, read_I2, read_R1, read_R2 string
	var strand_R1, strand_R2 string
	var qual_R1, qual_R2 string
	var index_p7, index_i7, index_p5, index_i5 = "", "", "", ""
	var index_1, index_2 string
	var to_write_R1, to_write_R2 string

	lengthP7 := LENGTHDIC["p7"]
	lengthP5 := LENGTHDIC["p5"]
	lengthI5 := LENGTHDIC["i5"]
	lengthI7 := LENGTHDIC["i7"]

	isP7 := lengthP7 > 0
	isP5 := lengthP5 > 0
	isI7 := lengthI7 > 0
	isI5 := lengthI5 > 0

	count := 0

	mainloop:
	for scannerI1.Scan() {
		scannerI2.Scan()
		scannerR1.Scan()
		scannerR2.Scan()

		id_I1 = scannerI1.Text()
		id_I2 = scannerI2.Text()
		id_R1 = scannerR1.Text()
		id_R2 = scannerR2.Text()

		if id_I1[0] != byte('@') || id_I2[0] != byte('@') ||
			id_R1[0] != byte('@') || id_R2[0] != byte('@') {
			continue mainloop
		}

		scannerI1.Scan()
		scannerI2.Scan()
		scannerR1.Scan()
		scannerR2.Scan()

		read_I1 = scannerI1.Text()
		read_I2 = scannerI2.Text()
		read_R1 = scannerR1.Text()
		read_R2 = scannerR2.Text()

		scannerI1.Scan()
		scannerI2.Scan()
		scannerR1.Scan()
		scannerR2.Scan()

		strand_R1 = scannerR1.Text()
		strand_R2 = scannerR2.Text()

		scannerI1.Scan()
		scannerI2.Scan()
		scannerR1.Scan()
		scannerR2.Scan()

		qual_R1 = scannerR1.Text()
		qual_R2 = scannerR2.Text()

		if TAGLENGTH > len(read_I1) || TAGLENGTH > len(read_I2) {
			fmt.Printf("#### error at read nb: %d\n", count)
			fmt.Printf("## error! ID I1: %s\n", id_I1)
			fmt.Printf("## error! ID I2: %s\n", id_I2)
			fmt.Printf("## error! ID I1: %s\n", id_R1)
			fmt.Printf("## error! ID I2: %s\n", id_R2)
			fmt.Printf("## error! read I1: %s\n", read_I1)
			fmt.Printf("## error! read I2: %s\n", read_I2)
			fmt.Printf("## error! read R1: %s\n", read_R1)
			fmt.Printf("## error! read R2: %s\n", read_R2)
			fmt.Printf("## error! strand R1: %s\n", strand_R1)
			fmt.Printf("## error! read I2: %s\n", strand_R2)

			if ERRORHANDLING != "" {
				switch ERRORHANDLING {

				case "return":
					fmt.Printf("error handling strategy: return at this stage and continue the pipeline\n")
					goto endfunc
				case "raise":
					err := fmt.Sprintf("#### error handling strategy: raise error #### \n")
					err += fmt.Sprintf("#### error at read nb: %d\n", count)
					if TAGLENGTH > len(read_I1) {
						err += fmt.Sprintf("error with id_I1: %s read:%s\n",
							id_I1, read_I1)
					} else {
						err += fmt.Sprintf("error with id_I2: %s read:%s\n",
							id_I2, read_I2)
					}
					panic(err)
				}
			}
		}

		if isP7 {
			index_p7 = read_I1[:lengthP7]
		}
		if isI7 {
			index_i7 = read_I1[len(read_I1)-lengthI7:]
		}
		if isI5 {
			index_i5 = read_I2[:lengthI5]
		}
		if isP5 {
			index_p5 = read_I2[len(read_I2)-lengthP5-SHIFT_P5:len(read_I2)-SHIFT_P5]
		}

		isValid, Replicate := checkIndexes(&index_p7, &index_i7, &index_p5, &index_i5)

		index = fmt.Sprintf("%s%s%s%s", index_p7, index_i7, index_i5, index_p5)

		if WRITE_LOGS {
			switch {
			case isValid:
				if _, isInside := logs[fmt.Sprintf("success_repl%d", Replicate)].dict[index]; !isInside {
					logs["stats"].dict[fmt.Sprintf("Number of cells repl. %d", Replicate)]++

					logsIndexCell[fmt.Sprintf("success_p5_repl%d", Replicate)].dict[index_p5]++
					logsIndexCell[fmt.Sprintf("success_p7_repl%d", Replicate)].dict[index_p7]++
					logsIndexCell[fmt.Sprintf("success_i5_repl%d", Replicate)].dict[index_i5]++
					logsIndexCell[fmt.Sprintf("success_i7_repl%d", Replicate)].dict[index_i7]++

				}
				logs[fmt.Sprintf("success_repl%d", Replicate)].dict[index]++
				logs["stats"].dict[fmt.Sprintf("Number of reads repl. %d", Replicate)]++
				logsIndexRead[fmt.Sprintf("success_p5_repl%d", Replicate)].dict[index_p5]++
				logsIndexRead[fmt.Sprintf("success_p7_repl%d", Replicate)].dict[index_p7]++
				logsIndexRead[fmt.Sprintf("success_i5_repl%d", Replicate)].dict[index_i5]++
				logsIndexRead[fmt.Sprintf("success_i7_repl%d", Replicate)].dict[index_i7]++

			default:
				if _, isInside := logs["fail"].dict[index]; !isInside {
					logs["stats"].dict["number of cells (FAIL)"]++

					logsIndexCell["fail_p5"].dict[index_p5]++
					logsIndexCell["fail_p7"].dict[index_p7]++
					logsIndexCell["fail_i5"].dict[index_i5]++
					logsIndexCell["fail_i7"].dict[index_i7]++
				}

				logs["fail"].dict[index]++
				logs["stats"].dict["number of reads (FAIL) "]++

				logsIndexRead["fail_p5"].dict[index_p5]++
				logsIndexRead["fail_p7"].dict[index_p7]++
				logsIndexRead["fail_i5"].dict[index_i5]++
				logsIndexRead["fail_i7"].dict[index_i7]++
			}

		}

		index_1 = fmt.Sprintf("@%s:%s", index, id_I1[1:])
		index_2 = fmt.Sprintf("@%s:%s", index, id_I2[1:])

		if !isValid {
			goto endmainloop
		}

		to_write_R1 = fmt.Sprintf("%s\n%s\n%s\n%s\n", index_1, read_R1, strand_R1, qual_R1)
		to_write_R2 = fmt.Sprintf("%s\n%s\n%s\n%s\n", index_2, read_R2, strand_R2, qual_R2)

		writerR1Dict[writerR1DictName[Replicate]].Write([]byte(to_write_R1))
		writerR2Dict[writerR2DictName[Replicate]].Write([]byte(to_write_R2))

		endmainloop:
		count++
		if maxNbReads != 0 && count >= maxNbReads {
			break mainloop
		}
	}
endfunc:
	for key, value := range logs {
		LOG_CHAN[key] <- *value
	}

	for key, value := range logsIndexCell {
		LOG_INDEX_CELL_CHAN[key] <- *value
	}

	for key, value := range logsIndexRead {
		LOG_INDEX_READ_CHAN[key] <- *value
	}
}


/*Check ... */
func Check(err error) {
	if err != nil {
		panic(err)
	}
}
