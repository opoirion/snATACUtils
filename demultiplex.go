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


import (
	utils "ATACdemultiplex/ATACdemultiplexUtils"
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




type StatsDict struct {
	dict map[string]int
	name string
}

func (statDict *StatsDict) Init() {
	statDict.dict = make(map[string]int)
}


func initChan(log_chan * map[string]chan StatsDict, log_type []string) {

	(*log_chan) = make(map[string] chan StatsDict)

	for _, value := range log_type {
		(*log_chan)[value] = make(chan StatsDict, NB_THREADS)
	}


}


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

	fmt.Printf("fastq file index 1 analyzed: %s\n", FASTQ_I1)
	fmt.Printf("fastq file index 2 analyzed: %s\n", FASTQ_I2)
	fmt.Printf("fastq file read file 1 analyzed: %s\n", FASTQ_R1)
	fmt.Printf("fastq file read file 2 analyzed: %s\n", FASTQ_R2)

	t_start := time.Now()

	switch {
	case NB_THREADS == 1:
		var waiting sync.WaitGroup
		waiting.Add(1)
		launchAnalysisOneFile(0, MAX_NB_READS, "", "", &waiting)

	case NB_THREADS < 1:
		log.Fatal(errors.New("threads should be >= 1!"))

	case NB_THREADS > 1:
		launchAnalysisMultipleFile("")
	}

	writeReport()

	t_diff := time.Now().Sub(t_start)
	fmt.Printf("demultiplexing finished in %f s\n", t_diff.Seconds())
}

func writeReport() {

	if !WRITE_LOGS{
		return
	}

	writeReportFromMultipleDict(&LOG_INDEX_CELL_CHAN, "index_cell")
	writeReportFromMultipleDict(&LOG_INDEX_READ_CHAN, "index_read")

	for key, log_chan := range LOG_CHAN {
		logs := extractDictFromChan(log_chan)
		filename := fmt.Sprintf("report%s_%s.log",
			OUTPUT_TAG_NAME, key)

		file, err := os.Create(filename)
		defer file.Close()
		Check(err)

		file.WriteString("#<key>\t<value>\n")

		for key, value := range logs {
			file.WriteString(fmt.Sprintf("%s\t%d\n", key, value))
		}

	}
}

func writeReportFromMultipleDict(channel * map[string]chan StatsDict, fname string) {
	filename := fmt.Sprintf("report%s_%s.log",
		OUTPUT_TAG_NAME, fname)
	file, err := os.Create(filename)

	defer file.Close()
	Check(err)


	filename = fmt.Sprintf("report%s_%s_fail.log",
		OUTPUT_TAG_NAME, fname)
	file_fail, err := os.Create(filename)

	defer file_fail.Close()
	Check(err)

	var fp *os.File

	for dict_type, log_chan := range *channel {

		switch{
		case strings.Contains(dict_type, "fail"):
			fp = file_fail
		default:
			fp = file
		}

		logs := extractDictFromChan(log_chan)
		fp.WriteString(fmt.Sprintf("#### %s\n", dict_type))

		for key, value := range logs {
			fp.WriteString(fmt.Sprintf("%s\t%s\t%d\n", dict_type, key, value))
		}
		fp.WriteString("\n")
	}
}

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

func launchAnalysisMultipleFile(path string) {

	var nb_reads int
	var chunk int

	switch {
	case MAX_NB_READS != 0:
		chunk = (MAX_NB_READS / NB_THREADS)
	default:
		fmt.Printf("computing number of lines...")
		nb_reads = countLine(FASTQ_I1, COMPRESSION_MODE) / 4
		fmt.Printf("estimated number of lines:%d\n", nb_reads)

		chunk = ((nb_reads / NB_THREADS) - (nb_reads/NB_THREADS)%4)
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
			"",
			fmt.Sprintf("index_%d.", index),
			&waiting)
		startingRead += chunk
		index += 1
	}

	waiting.Wait()

	output_R1 := fmt.Sprintf("%s%s%s.demultiplexed.repl1.bz2", path,
		strings.TrimSuffix(FASTQ_R1, ".bz2"),
		OUTPUT_TAG_NAME)
	output_R2 := fmt.Sprintf("%s%s%s.demultiplexed.repl1.bz2", path,
		strings.TrimSuffix(FASTQ_R2, ".bz2"),
		OUTPUT_TAG_NAME)

	cmd_1 := fmt.Sprintf("cat index_*demultiplexed.repl1.bz2 > %s", output_R1)
	cmd_2 := fmt.Sprintf("cat index_*demultiplexed.repl1.bz2 > %s", output_R2)

	fmt.Printf("concatenating read 1 files...\n")
	fmt.Printf("%s\n", cmd_1)

	utils.ExceCmd(cmd_1)
	fmt.Printf("concatenating read 2 files...\n")
	utils.ExceCmd(cmd_2)

	if INDEX_REPLICATE_R2 != "" {
		output_R1 := fmt.Sprintf("%s%s%s.demultiplexed.repl2.bz2", path,
			strings.TrimSuffix(FASTQ_R1, ".bz2"),
			OUTPUT_TAG_NAME)
		output_R2 := fmt.Sprintf("%s%s%s.demultiplexed.repl2.bz2", path,
			strings.TrimSuffix(FASTQ_R2, ".bz2"),
			OUTPUT_TAG_NAME)

		cmd_1 := fmt.Sprintf("cat index_*demultiplexed.repl2.bz2 > %s", output_R1)
		cmd_2 := fmt.Sprintf("cat index_*demultiplexed.repl2.bz2 > %s", output_R2)

		fmt.Printf("concatenating read 1 files...\n")
		fmt.Printf("%s\n", cmd_1)

		utils.ExceCmd(cmd_1)
		fmt.Printf("concatenating read 2 files...\n")
		utils.ExceCmd(cmd_2)
	}

	cmd := fmt.Sprintf("rm index_*demultiplexed.repl*.bz2 ")

	fmt.Printf("removing read index files...\n")
	utils.ExceCmd(cmd)
}

func initLog(log_type []string) (logs map[string]*StatsDict) {

	logs = make(map[string]*StatsDict)

	for _, value := range log_type {
			logs[value] = &StatsDict{name:value}
			logs[value].Init()
	}

	return logs
}

func launchAnalysisOneFile(
	startingRead int,
	max_nb_reads int,
	path string,
	index string,
	waiting *sync.WaitGroup) {

	defer waiting.Done()

	logs := initLog(LOG_TYPE)
	logs_index_read := initLog(LOG_INDEX_TYPE)
	logs_index_cell := initLog(LOG_INDEX_TYPE)

	var scanner_I1 * bufio.Scanner
	var scanner_I2 * bufio.Scanner
	var scanner_R1 * bufio.Scanner
	var scanner_R2 * bufio.Scanner
	var bzip_R1_repl1  io.WriteCloser
	var bzip_R2_repl1  io.WriteCloser
	var bzip_R1_repl2  io.WriteCloser
	var bzip_R2_repl2  io.WriteCloser

	var file_I1 * os.File
	var file_I2 * os.File
	var file_R1 * os.File
	var file_R2 * os.File

	output_R1_repl1 := fmt.Sprintf("%s%s%s%s.demultiplexed.repl1.bz2",
		path, index, strings.TrimSuffix(FASTQ_R1, ".bz2"),
		OUTPUT_TAG_NAME)
	output_R2_repl1 := fmt.Sprintf("%s%s%s%s.demultiplexed.repl1.bz2", path, index,
		strings.TrimSuffix(FASTQ_R2, ".bz2"),
		OUTPUT_TAG_NAME)

	output_R1_repl2 := fmt.Sprintf("%s%s%s%s.demultiplexed.repl2.bz2",
		path, index, strings.TrimSuffix(FASTQ_R1, ".bz2"),
		OUTPUT_TAG_NAME)
	output_R2_repl2 := fmt.Sprintf("%s%s%s%s.demultiplexed.repl2.bz2", path, index,
		strings.TrimSuffix(FASTQ_R2, ".bz2"),
		OUTPUT_TAG_NAME)

	switch  {
	case USE_BZIP_GO_LIBRARY:
		scanner_I1, file_I1 = utils.ReturnReaderForBzipfilePureGo(FASTQ_I1, startingRead * 4)
		scanner_I2, file_I2 = utils.ReturnReaderForBzipfilePureGo(FASTQ_I2, startingRead * 4)
		scanner_R1, file_R1 = utils.ReturnReaderForBzipfilePureGo(FASTQ_R1, startingRead * 4)
		scanner_R2, file_R2 = utils.ReturnReaderForBzipfilePureGo(FASTQ_R2, startingRead * 4)

		GUESS_NB_LINES = false

		bzip_R1_repl1 = utils.ReturnWriterForBzipfilePureGo(output_R1_repl1)
		bzip_R2_repl1 = utils.ReturnWriterForBzipfilePureGo(output_R2_repl1)

		if INDEX_REPLICATE_R2 != "" {
			bzip_R1_repl2 = utils.ReturnWriterForBzipfilePureGo(output_R1_repl2)
			bzip_R2_repl2 = utils.ReturnWriterForBzipfilePureGo(output_R2_repl2)
			defer bzip_R1_repl2.Close()
			defer bzip_R2_repl2.Close()
		}

	default:
		scanner_I1, file_I1 = utils.ReturnReaderForBzipfile(FASTQ_I1, startingRead * 4)
		scanner_I2, file_I2 = utils.ReturnReaderForBzipfile(FASTQ_I2, startingRead * 4)
		scanner_R1, file_R1 = utils.ReturnReaderForBzipfile(FASTQ_R1, startingRead * 4)
		scanner_R2, file_R2 = utils.ReturnReaderForBzipfile(FASTQ_R2, startingRead * 4)
		bzip_R1_repl1 = utils.ReturnWriterForBzipfile(output_R1_repl1, COMPRESSION_MODE)
		bzip_R2_repl1 = utils.ReturnWriterForBzipfile(output_R2_repl1, COMPRESSION_MODE)

		if INDEX_REPLICATE_R2 != "" {
			bzip_R1_repl2 = utils.ReturnWriterForBzipfile(output_R1_repl2, COMPRESSION_MODE)
			bzip_R2_repl2 = utils.ReturnWriterForBzipfile(output_R2_repl2, COMPRESSION_MODE)
			defer bzip_R1_repl2.Close()
			defer bzip_R2_repl2.Close()
		}

	}

	defer file_I1.Close()
	defer file_I2.Close()
	defer file_R1.Close()
	defer file_R2.Close()
	defer bzip_R1_repl1.Close()
	defer bzip_R2_repl1.Close()

	var id_I1, id_I2, id_R1, id_R2 string
	var read_I1, read_I2, read_R1, read_R2 string
	var strand_R1, strand_R2 string
	var qual_R1, qual_R2 string
	var index_p7, index_i7, index_p5, index_i5 string
	var index_1, index_2 string
	var to_write_R1, to_write_R2 string

	count := 0

	mainloop:
	for scanner_I1.Scan() {
		scanner_I2.Scan()
		scanner_R1.Scan()
		scanner_R2.Scan()

		id_I1 = scanner_I1.Text()
		id_I2 = scanner_I2.Text()
		id_R1 = scanner_R1.Text()
		id_R2 = scanner_R2.Text()

		if id_I1[0] != byte('@') || id_I2[0] != byte('@') ||
			id_R1[0] != byte('@') || id_R2[0] != byte('@') {
			continue mainloop
		}

		scanner_I1.Scan()
		scanner_I2.Scan()
		scanner_R1.Scan()
		scanner_R2.Scan()

		read_I1 = scanner_I1.Text()
		read_I2 = scanner_I2.Text()
		read_R1 = scanner_R1.Text()
		read_R2 = scanner_R2.Text()

		scanner_I1.Scan()
		scanner_I2.Scan()
		scanner_R1.Scan()
		scanner_R2.Scan()

		strand_R1 = scanner_R1.Text()
		strand_R2 = scanner_R2.Text()

		scanner_I1.Scan()
		scanner_I2.Scan()
		scanner_R1.Scan()
		scanner_R2.Scan()

		qual_R1 = scanner_R1.Text()
		qual_R2 = scanner_R2.Text()

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

					logs_index_cell["success_p5_repl1"].dict[index_p5]++
					logs_index_cell["success_p7_repl1"].dict[index_p7]++
					logs_index_cell["success_i5_repl1"].dict[index_i5]++
					logs_index_cell["success_i7_repl1"].dict[index_i7]++

				}
				logs["success_repl_1"].dict[index] += 1
				logs["stats"].dict["Number of reads repl. 1"]++

				logs_index_read["success_p5_repl1"].dict[index_p5]++
				logs_index_read["success_p7_repl1"].dict[index_p7]++
				logs_index_read["success_i5_repl1"].dict[index_i5]++
				logs_index_read["success_i7_repl1"].dict[index_i7]++

			case isValid && Replicate == 2:
				if _, isInside := logs["success_repl_2"].dict[index]; !isInside {
					logs["stats"].dict["Number of cells repl. 2"]++

					logs_index_cell["success_p5_repl2"].dict[index_p5]++
					logs_index_cell["success_p7_repl2"].dict[index_p7]++
					logs_index_cell["success_i5_repl2"].dict[index_i5]++
					logs_index_cell["success_i7_repl2"].dict[index_i7]++
				}

				logs["success_repl_2"].dict[index] += 1
				logs["stats"].dict["Number of reads repl. 2"]++

				logs_index_read["success_p5_repl2"].dict[index_p5]++
				logs_index_read["success_p7_repl2"].dict[index_p7]++
				logs_index_read["success_i5_repl2"].dict[index_i5]++
				logs_index_read["success_i7_repl2"].dict[index_i7]++

			default:
				if _, isInside := logs["fail"].dict[index]; !isInside {
					logs["stats"].dict["number of cells (FAIL)"]++

					logs_index_cell["fail_p5"].dict[index_p5]++
					logs_index_cell["fail_p7"].dict[index_p7]++
					logs_index_cell["fail_i5"].dict[index_i5]++
					logs_index_cell["fail_i7"].dict[index_i7]++
				}

				logs["fail"].dict[index] += 1
				logs["stats"].dict["number of reads (FAIL) "]++

				logs_index_read["fail_p5"].dict[index_p5]++
				logs_index_read["fail_p7"].dict[index_p7]++
				logs_index_read["fail_i5"].dict[index_i5]++
				logs_index_read["fail_i7"].dict[index_i7]++
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
			bzip_R1_repl1.Write([]byte(to_write_R1))
			bzip_R2_repl1.Write([]byte(to_write_R2))
		case Replicate == 2:
			bzip_R1_repl2.Write([]byte(to_write_R1))
			bzip_R2_repl2.Write([]byte(to_write_R2))
		default:
			log.Fatal("Error wrong replicate!")
		}

		endmainloop:
		count += 1
		if max_nb_reads != 0 && count >= max_nb_reads {
			break mainloop
		}
	}

	for key, value := range logs {
		LOG_CHAN[key] <- *value
	}

	for key, value := range logs_index_cell {
		LOG_INDEX_CELL_CHAN[key] <- *value
	}

	for key, value := range logs_index_read {
		LOG_INDEX_READ_CHAN[key] <- *value
	}

}


func Check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
