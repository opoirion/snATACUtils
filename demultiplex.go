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

var FASTQ_R1 string
var FASTQ_R2 string
var FASTQ_I1 string
var FASTQ_I2 string
var NB_THREADS int
var TAGLENGTH int
var MAX_NB_READS int
var COMPRESSION_MODE int
var USE_BZIP_GO_LIBRARY bool
var GUESS_NB_LINES bool
var WRITE_LOGS bool
var OUTPUT_TAG_NAME string
var INDEX_REPLICATE_R1 string
var INDEX_REPLICATE_R2 string
var INDEX_NO_REPLICATE string
var OUTPUT_PATH string
var MAX_NB_MISTAKE int

var INDEX_R1_DICT map[string]map[string]bool
var INDEX_R2_DICT map[string]map[string]bool
var INDEX_NO_DICT map[string]map[string]bool


var LOG_CHAN map[string]chan StatsDict
var LOG_INDEX_READ_CHAN map[string]chan StatsDict
var LOG_INDEX_CELL_CHAN map[string]chan StatsDict

var LOG_TYPE  = []string {
	"stats",
	"success_repl_1",
	"success_repl_2",
	"fail",
}

var LOG_INDEX_TYPE  = []string {
	"success_p5_repl1",
	"success_p7_repl1",
	"success_i5_repl1",
	"success_i7_repl1",

	"success_p5_repl2",
	"success_p7_repl2",
	"success_i5_repl2",
	"success_i7_repl2",

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
	flag.IntVar(&MAX_NB_MISTAKE, "max_nb_mistake", 2, "Maximum number of mistakes allowed to assign a reference read id (default 2)")
	flag.StringVar(&OUTPUT_TAG_NAME, "output_tag_name", "", "tag for the output file names (default None)")
	flag.BoolVar(&USE_BZIP_GO_LIBRARY, "use_bzip2_go_lib", false, "use bzip2 go library instead of native C lib (slower)")
	flag.BoolVar(&WRITE_LOGS, "write_logs", false, "write logs (might slower the execution time)")
	flag.BoolVar(&GUESS_NB_LINES, "guess_nb_lines", false, "guess automatically position of the lines (for mulithread). May be not safe in some situation")

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

	flag.Parse()

	if INDEX_REPLICATE_R1 != "" || INDEX_REPLICATE_R2 != "" {
		if  INDEX_NO_REPLICATE != "" {
			log.Fatal("Cannot set up index_no_replicate with either index_replicate_r1 or index_replicate_r2")
		}
		if INDEX_REPLICATE_R1 == "" || INDEX_REPLICATE_R2 == "" {
			log.Fatal("Both index_replicate_r1 and index_replicate_r2 should be set up!")
		}
	}

	initChan(&LOG_CHAN, LOG_TYPE)
	initChan(&LOG_INDEX_READ_CHAN, LOG_INDEX_TYPE)
	initChan(&LOG_INDEX_CELL_CHAN, LOG_INDEX_TYPE)

	loadIndexes(INDEX_NO_REPLICATE, &INDEX_NO_DICT)
	loadIndexes(INDEX_REPLICATE_R1, &INDEX_R1_DICT)
	loadIndexes(INDEX_REPLICATE_R2, &INDEX_R2_DICT)

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

		for key, value := range logs {
			file.WriteString(fmt.Sprintf("%s\t%d\n", key, value))
		}

	}
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
		default:
			fp = file
		}

		logs := extractDictFromChan(logChan)
		fp.WriteString(fmt.Sprintf("#### %s\n", dictType))

		for key, value := range logs {
			fp.WriteString(fmt.Sprintf("%s\t%s\t%d\n", dictType, key, value))
		}
		fp.WriteString("\n")
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
		fmt.Printf("computing number of lines...")
		nbReads = countLine(FASTQ_I1, COMPRESSION_MODE) / 4
		fmt.Printf("estimated number of lines:%d\n", nbReads)

		chunk = ((nbReads / NB_THREADS) - (nbReads/NB_THREADS)%4)
	}

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

	outputR1 := fmt.Sprintf("%s%s%s.demultiplexed.repl1.fastq.%s", OUTPUT_PATH,
		strings.TrimSuffix(filenameR1, fmt.Sprintf(".fastq.%s", ext)),
		OUTPUT_TAG_NAME, ext)
	outputR2 := fmt.Sprintf("%s%s%s.demultiplexed.repl1.fastq.%s", OUTPUT_PATH,
		strings.TrimSuffix(filenameR2, fmt.Sprintf(".fastq.%s", ext)),
		OUTPUT_TAG_NAME, ext)

	cmd1 := fmt.Sprintf("cat %sindex_*demultiplexed.repl1.fastq.bz2 > %s",
		OUTPUT_PATH, outputR1)
	cmd2 := fmt.Sprintf("cat %sindex_*demultiplexed.repl1.fastq.bz2 > %s",
		OUTPUT_PATH, outputR2)

	fmt.Printf("concatenating repl. 1 read 1 files...\n")
	utils.ExceCmd(cmd1)
	fmt.Printf("concatenating repl. 1 read 2 files...\n")
	utils.ExceCmd(cmd2)

	if INDEX_REPLICATE_R2 != "" {
		outputR1 := fmt.Sprintf("%s%s%s.demultiplexed.repl2.fastq.%s", OUTPUT_PATH,
			strings.TrimSuffix(filenameR1, fmt.Sprintf(".fastq.%s", ext)),
			OUTPUT_TAG_NAME, ext)
		outputR2 := fmt.Sprintf("%s%s%s.demultiplexed.repl2.fastq.%s", OUTPUT_PATH,
			strings.TrimSuffix(filenameR2, fmt.Sprintf(".fastq.%s", ext)),
			OUTPUT_TAG_NAME, ext)

		cmd1 := fmt.Sprintf("cat %sindex_*demultiplexed.repl2.fastq.bz2 > %s",
			OUTPUT_PATH, outputR1)
		cmd2 := fmt.Sprintf("cat %sindex_*demultiplexed.repl2.fastq.bz2 > %s",
			OUTPUT_PATH, outputR2)

		fmt.Printf("concatenating repl. 2 read 1 files...\n")
		fmt.Printf("%s\n", cmd1)

		utils.ExceCmd(cmd1)
		fmt.Printf("concatenating repl. 2 read 2 files...\n")
		utils.ExceCmd(cmd2)
	}

	cmd := fmt.Sprintf("rm %sindex_*demultiplexed.repl*.bz2 ", OUTPUT_PATH)

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
	var bzipR1repl1  io.WriteCloser
	var bzipR2repl1  io.WriteCloser
	var bzipR1repl2  io.WriteCloser
	var bzipR2repl2  io.WriteCloser

	var fileI1 * os.File
	var fileI2 * os.File
	var fileR1 * os.File
	var fileR2 * os.File

	_, filenameR1 := pathutils.Split(FASTQ_R1)
	_, filenameR2 := pathutils.Split(FASTQ_R2)

	ext := path.Ext(filenameR1)

	outputR1Repl1 := fmt.Sprintf("%s%s%s%s.demultiplexed.repl1.fastq.%s",
		OUTPUT_PATH, index,
		strings.TrimSuffix(filenameR1, fmt.Sprintf(".fastq.%s", ext)),
		OUTPUT_TAG_NAME, ext)
	outputR2Repl1 := fmt.Sprintf("%s%s%s%s.demultiplexed.repl1.fastq.%s", OUTPUT_PATH, index,
		strings.TrimSuffix(filenameR2, fmt.Sprintf(".fastq.%s", ext)),
		OUTPUT_TAG_NAME, ext)

	outputR1Repl2 := fmt.Sprintf("%s%s%s%s.demultiplexed.repl2.fastq.%s",
		OUTPUT_PATH, index,
		strings.TrimSuffix(filenameR1, fmt.Sprintf(".fastq.%s", ext)),
		OUTPUT_TAG_NAME, ext)
	outputR2Repl2 := fmt.Sprintf("%s%s%s%s.demultiplexed.repl2.fastq.%s", OUTPUT_PATH, index,
		strings.TrimSuffix(filenameR2, fmt.Sprintf(".fastq.%s", ext)),
		OUTPUT_TAG_NAME, ext)


	scannerI1, fileI1 = utils.ReturnReader(FASTQ_I1, startingRead * 4, USE_BZIP_GO_LIBRARY)
	scannerI2, fileI2 = utils.ReturnReader(FASTQ_I2, startingRead * 4, USE_BZIP_GO_LIBRARY)
	scannerR1, fileR1 = utils.ReturnReader(FASTQ_R1, startingRead * 4, USE_BZIP_GO_LIBRARY)
	scannerR2, fileR2 = utils.ReturnReader(FASTQ_R2, startingRead * 4, USE_BZIP_GO_LIBRARY)

	GUESS_NB_LINES = false

	bzipR1repl1 = utils.ReturnWriter(outputR1Repl1, COMPRESSION_MODE,
		USE_BZIP_GO_LIBRARY)
	bzipR2repl1 = utils.ReturnWriter(outputR2Repl1, COMPRESSION_MODE,
		USE_BZIP_GO_LIBRARY)

	if INDEX_REPLICATE_R2 != "" {
		bzipR1repl2 = utils.ReturnWriter(outputR1Repl2, COMPRESSION_MODE,
			USE_BZIP_GO_LIBRARY)
		bzipR2repl2 = utils.ReturnWriter(outputR2Repl2, COMPRESSION_MODE,
			USE_BZIP_GO_LIBRARY)
		defer bzipR1repl2.Close()
		defer bzipR2repl2.Close()
	}

	defer fileI1.Close()
	defer fileI2.Close()
	defer fileR1.Close()
	defer fileR2.Close()
	defer bzipR1repl1.Close()
	defer bzipR2repl1.Close()

	var id_I1, id_I2, id_R1, id_R2 string
	var read_I1, read_I2, read_R1, read_R2 string
	var strand_R1, strand_R2 string
	var qual_R1, qual_R2 string
	var index_p7, index_i7, index_p5, index_i5 string
	var index_1, index_2 string
	var to_write_R1, to_write_R2 string

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

		index_p7 = read_I1[:TAGLENGTH]
		index_i7 = read_I1[len(read_I1)-TAGLENGTH:]
		index_i5 = read_I2[:TAGLENGTH]
		index_p5 = read_I2[len(read_I2)-TAGLENGTH:]

		isValid, Replicate := checkIndexes(&index_p7, &index_i7, &index_p5, &index_i5)

		index = fmt.Sprintf("%s%s%s%s", index_p7, index_i7, index_p5, index_i5)

		if WRITE_LOGS {
			switch {
			case isValid && Replicate == 1:
				if _, isInside := logs["success_repl_1"].dict[index]; !isInside {
					logs["stats"].dict["Number of cells repl. 1"]++

					logsIndexCell["success_p5_repl1"].dict[index_p5]++
					logsIndexCell["success_p7_repl1"].dict[index_p7]++
					logsIndexCell["success_i5_repl1"].dict[index_i5]++
					logsIndexCell["success_i7_repl1"].dict[index_i7]++

				}
				logs["success_repl_1"].dict[index]++
				logs["stats"].dict["Number of reads repl. 1"]++

				logsIndexRead["success_p5_repl1"].dict[index_p5]++
				logsIndexRead["success_p7_repl1"].dict[index_p7]++
				logsIndexRead["success_i5_repl1"].dict[index_i5]++
				logsIndexRead["success_i7_repl1"].dict[index_i7]++

			case isValid && Replicate == 2:
				if _, isInside := logs["success_repl_2"].dict[index]; !isInside {
					logs["stats"].dict["Number of cells repl. 2"]++

					logsIndexCell["success_p5_repl2"].dict[index_p5]++
					logsIndexCell["success_p7_repl2"].dict[index_p7]++
					logsIndexCell["success_i5_repl2"].dict[index_i5]++
					logsIndexCell["success_i7_repl2"].dict[index_i7]++
				}

				logs["success_repl_2"].dict[index]++
				logs["stats"].dict["Number of reads repl. 2"]++

				logsIndexRead["success_p5_repl2"].dict[index_p5]++
				logsIndexRead["success_p7_repl2"].dict[index_p7]++
				logsIndexRead["success_i5_repl2"].dict[index_i5]++
				logsIndexRead["success_i7_repl2"].dict[index_i7]++

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

		switch {
		case Replicate == 1:
			bzipR1repl1.Write([]byte(to_write_R1))
			bzipR2repl1.Write([]byte(to_write_R2))
		case Replicate == 2:
			bzipR1repl2.Write([]byte(to_write_R1))
			bzipR2repl2.Write([]byte(to_write_R2))
		default:
			log.Fatal("Error wrong replicate!")
		}

		endmainloop:
		count++
		if maxNbReads != 0 && count >= maxNbReads {
			break mainloop
		}
	}

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


/* */
func Check(err error) {
	if err != nil {
		panic(err)
	}
}
