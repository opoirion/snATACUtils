package main

import "fmt"
import "flag"
import "os"
import "log"
import "bufio"
import "strings"
import "github.com/dsnet/compress/bzip2"
import "time"
import "errors"
import "os/exec"
// import "sync"

var FASTQ_R1 string
var FASTQ_R2 string
var FASTQ_I1 string
var FASTQ_I2 string
var NB_THREADS int
var NB_THREADS_CHUNK int
var TAGLENGTH int
var MAX_NB_READS int

type ReaderStatus struct {
	isEOF bool
}


func main() {
	// counter := make(map[string]int)
	flag.StringVar(&FASTQ_I1, "fastq_I1", "", "fastq index file index paired read 1")
	flag.StringVar(&FASTQ_I2, "fastq_I2", "", "fastq index file index paired read 2")
	flag.StringVar(&FASTQ_R1, "fastq_R1", "", "fastq read file index paired read 1")
	flag.StringVar(&FASTQ_R2, "fastq_R2", "", "fastq read file index paired read 2")
	flag.IntVar(&NB_THREADS, "nbThreads", 1, "number of threads to use")
	flag.IntVar(&NB_THREADS_CHUNK, "nbThreadsChunk", 100000, `Chunck to be processed by each threads
(used only with NB_THREADS > 1)`)
	flag.IntVar(&TAGLENGTH, "taglength", 8,
		`<OPTIONAL> number of nucleotides to consider at the end
 and begining (default 8)`)
	flag.IntVar(&MAX_NB_READS, "max_nb_reads", 0,
		"<OPTIONAL> max number of reads to process (default 0 => None)")

	flag.Parse()

	fmt.Printf("fastq file index 1 analyzed: %s\n", FASTQ_I1)
	fmt.Printf("fastq file index 2 analyzed: %s\n", FASTQ_I2)
	fmt.Printf("fastq file read file 1 analyzed: %s\n", FASTQ_R1)
	fmt.Printf("fastq file read file 2 analyzed: %s\n", FASTQ_R2)

	t_start := time.Now().Second()

	switch {
	case NB_THREADS == 1:
		readerStatusChan := make(chan *ReaderStatus, 1)
		launchAnalysisOneFile(0, MAX_NB_READS, "", "",
			&ReaderStatus{isEOF:false}, readerStatusChan)
	case NB_THREADS < 1:
		log.Fatal(errors.New("threads should be >= 1!"))

	case NB_THREADS > 1:
		launchAnalysisMultipleFile("")
	}

	t_end := time.Now().Second()
	fmt.Printf("demultiplexing finished in %d s\n", t_end - t_start)
}


func launchAnalysisMultipleFile(path string) {
	readerStatusChan := make(chan *ReaderStatus, NB_THREADS)

	for i := 0; i < NB_THREADS; i++ {
		readerStatusChan <- &ReaderStatus{isEOF:false}
	}


	startingLine := 0
	index := 0
	nb_threads_chunk := NB_THREADS_CHUNK - NB_THREADS_CHUNK%4
	nbFinished := 0

	loop:
	for {
		select {
		case readerStatus := <-readerStatusChan:
			switch {
			case  MAX_NB_READS != 0 && startingLine >= MAX_NB_READS:
				nbFinished++
			case readerStatus.isEOF == false:
				go launchAnalysisOneFile(
					startingLine, nb_threads_chunk,
					"",
					fmt.Sprintf("index_%d", index),
					readerStatus,
					readerStatusChan,
					)

				index++
				startingLine += nb_threads_chunk

			case readerStatus.isEOF == true:
				fmt.Printf("One finished!\n")
				nbFinished++
				fmt.Printf("%d chunk read done\n", nb_threads_chunk)
			}

		default:
			switch {
			case nbFinished < NB_THREADS:
				time.Sleep(500 * time.Millisecond)
			default:
				fmt.Printf("break\n")
				break loop
			}

		}
	}

	// output_R1 := fmt.Sprintf("%s%s.demultiplexed.R1.bz2", path, strings.TrimSuffix(FASTQ_R1, ".bz2"))
	// // output_R2 := fmt.Sprintf("%s%s.demultiplexed.R2.bz2", path, strings.TrimSuffix(FASTQ_R2, ".bz2"))

	// cmd := fmt.Sprintf("cat index_*demultiplexed*R1.bz2 > %s", output_R1)
	// fmt.Printf("%s\n", cmd)

	// exe_cmd(cmd)
}


func exe_cmd(cmd string) {
  fmt.Println("command is ",cmd)
  // splitting head => g++ parts => rest of the command
  parts := strings.Fields(cmd)
  head := parts[0]
  parts = parts[1:len(parts)]

  err := exec.Command(head, parts...).Run()
  if err != nil {
    fmt.Printf("%s", err)
  }
}

func launchAnalysisOneFile(
	startingLine int,
	max_nb_reads int,
	path string,
	index string,
	readerStatus *ReaderStatus,
	readerStatusChan chan *ReaderStatus)  {

	if startingLine%4 != 0 {
		log.Fatal(errors.New("startingLine is not a multiple of 4!"))
	}

	//opening index 1
	scanner_I1 := return_reader_for_bzipfile(FASTQ_I1)

	//opening index 2
	scanner_I2 := return_reader_for_bzipfile(FASTQ_I2)

	//opening fastq read file 1
	scanner_R1 := return_reader_for_bzipfile(FASTQ_R1)

	//opening fastq read file 2
	scanner_R2 := return_reader_for_bzipfile(FASTQ_R2)

	output_R1 := fmt.Sprintf("%s%s.demultiplexed%s.R1.bz2", path, index, strings.TrimSuffix(FASTQ_R1, ".bz2"))
	output_R2 := fmt.Sprintf("%s%s.demultiplexed%s.R2.bz2", path, index, strings.TrimSuffix(FASTQ_R2, ".bz2"))

	bzip_R1 := return_writer_for_bzipfile(output_R1)
	bzip_R2 := return_writer_for_bzipfile(output_R2)

	defer bzip_R1.Close()
	defer bzip_R2.Close()

	count := 0
	for lineCount := 0; lineCount < startingLine; lineCount++ {
		statusI1 := scanner_I1.Scan()
		statusI2 := scanner_I2.Scan()
		statusR1 := scanner_R1.Scan()
		statusR2 := scanner_R2.Scan()

		if statusI1 == false || statusI2 == false || statusR1 == false || statusR2 == false {
			readerStatus.isEOF = true
			readerStatusChan <- readerStatus
			return
		}
	}

	interupted := false

	for scanner_I1.Scan() {
		scanner_I2.Scan()
		scanner_R1.Scan()
		scanner_R2.Scan()

		id_I1 := scanner_I1.Text()
		id_I2 := scanner_I2.Text()

		scanner_I1.Scan()
		scanner_I2.Scan()
		scanner_R1.Scan()
		scanner_R2.Scan()

		read_I1 := scanner_I1.Text()
		read_I2 := scanner_I2.Text()
		read_R1 := scanner_R1.Text()
		read_R2 := scanner_R2.Text()

		scanner_I1.Scan()
		scanner_I2.Scan()
		scanner_R1.Scan()
		scanner_R2.Scan()

		strand_R1 := scanner_R1.Text()
		strand_R2 := scanner_R2.Text()

		scanner_I1.Scan()
		scanner_I2.Scan()
		scanner_R1.Scan()
		scanner_R2.Scan()

		qual_R1 := scanner_R1.Text()
		qual_R2 := scanner_R2.Text()

		index_p7 := read_I1[:TAGLENGTH]
		index_i7 := read_I1[len(read_I1)-TAGLENGTH:]
		index_p5 := read_I2[:TAGLENGTH]
		index_i5 := read_I2[len(read_I2)-TAGLENGTH:]

		index := fmt.Sprintf("%s%s%s%s", index_p7, index_i7, index_p5, index_i5)
		index_1 := fmt.Sprintf("@%s:%s", index, id_I1[1:])
		index_2 := fmt.Sprintf("@%s:%s", index, id_I2[1:])

		to_write_R1 := fmt.Sprintf("%s\n%s\n%s\n%s\n", index_1,read_R1, strand_R1, qual_R1)
		to_write_R2 := fmt.Sprintf("%s\n%s\n%s\n%s\n", index_2, read_R2, strand_R2, qual_R2)

		bzip_R1.Write([]byte(to_write_R1))
		bzip_R2.Write([]byte(to_write_R2))

		count += 1

		if max_nb_reads !=0 && count >= max_nb_reads {
			interupted = true
			break
		}
	}

	if interupted == false {
		readerStatus.isEOF = true
	}

	readerStatusChan <- readerStatus
	return
}

func check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}

func return_writer_for_bzipfile(fname string) *bzip2.Writer {
	output_file, err := os.Create(fname)
	check(err)
	bzip_file, err := bzip2.NewWriter(output_file, new(bzip2.WriterConfig))
	check(err)

	return bzip_file
}

func return_reader_for_bzipfile(fname string) *bufio.Scanner {
	file_open, err := os.OpenFile(fname, 0, 0)

	if err != nil {
		log.Fatal(err)
	}
	config := new(bzip2.ReaderConfig)

	reader_os := bufio.NewReader(file_open)
	reader_bzip, err := bzip2.NewReader(reader_os, config)
	check(err)
	reader_os2 := bufio.NewReader(reader_bzip)
	bzip_scanner := bufio.NewScanner(reader_os2)

	return bzip_scanner
}
