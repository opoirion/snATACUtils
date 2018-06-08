package main

import "fmt"
import "flag"
import "os"
import "log"
import "bufio"
import "strings"
import "github.com/dsnet/compress/bzip2"
import "time"

var FASTQ_R1 string
var FASTQ_R2 string
var FASTQ_I1 string
var FASTQ_I2 string
var TAGLENGTH int
var MAX_NB_READS int

func main() {
	// counter := make(map[string]int)
	flag.StringVar(&FASTQ_I1, "fastq_I1", "", "fastq index file index paired read 1")
	flag.StringVar(&FASTQ_I2, "fastq_I2", "", "fastq index file index paired read 2")
	flag.StringVar(&FASTQ_R1, "fastq_R1", "", "fastq read file index paired read 1")
	flag.StringVar(&FASTQ_R2, "fastq_R2", "", "fastq read file index paired read 2")
	flag.IntVar(&TAGLENGTH, "taglength", 8, "<OPTIONAL> number of nucleotides to consider at the end and begining (default 8)")
	flag.IntVar(&MAX_NB_READS, "max_nb_reads", 0, "<OPTIONAL> max number of reads to process (default 0 => None)")

	flag.Parse()

	t_start := time.Now().Second()

	//opening index 1
	fmt.Printf("fastq file index 1 analyzed: %s\n", FASTQ_I1)
	scanner_I1 := return_reader_for_bzipfile(FASTQ_I1)

	//opening index 2
	fmt.Printf("fastq file index 2 analyzed: %s\n", FASTQ_I2)
	scanner_I2 := return_reader_for_bzipfile(FASTQ_I2)

	//opening fastq read file 1
	fmt.Printf("fastq file index 2 analyzed: %s\n", FASTQ_R1)
	scanner_R1 := return_reader_for_bzipfile(FASTQ_R1)

	//opening fastq read file 2
	fmt.Printf("fastq file index 2 analyzed: %s\n", FASTQ_R2)
	scanner_R2 := return_reader_for_bzipfile(FASTQ_R2)

	output_I1 := fmt.Sprintf("%s.demultiplexed.bz2", strings.SplitN(FASTQ_I1, ".", 2)[0])
	output_I2 := fmt.Sprintf("%s.demultiplexed.bz2", strings.SplitN(FASTQ_I2, ".", 2)[0])
	output_R1 := fmt.Sprintf("%s.demultiplexed.bz2", strings.SplitN(FASTQ_R1, ".", 2)[0])
	output_R2 := fmt.Sprintf("%s.demultiplexed.bz2", strings.SplitN(FASTQ_R2, ".", 2)[0])

	fmt.Printf("output file for I1: %s\n", output_I1)
	fmt.Printf("output file for R1: %s\n", output_R1)
	fmt.Printf("output file for I2: %s\n", output_I2)
	fmt.Printf("output file for I2: %s\n", output_R2)

	bzip_I1 := return_writer_for_bzipfile(output_I1)
	bzip_I2 := return_writer_for_bzipfile(output_I2)
	bzip_R1 := return_writer_for_bzipfile(output_R1)
	bzip_R2 := return_writer_for_bzipfile(output_R2)

	defer bzip_I1.Close()
	defer bzip_I2.Close()
	defer bzip_R1.Close()
	defer bzip_R2.Close()

	count := 0

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

		strand_I1 := scanner_I1.Text()
		strand_I2 := scanner_I2.Text()
		strand_R1 := scanner_R1.Text()
		strand_R2 := scanner_R2.Text()

		scanner_I1.Scan()
		scanner_I2.Scan()
		scanner_R1.Scan()
		scanner_R2.Scan()

		qual_I1 := scanner_I1.Text()
		qual_I2 := scanner_I2.Text()
		qual_R1 := scanner_R1.Text()
		qual_R2 := scanner_R2.Text()

		index_p7 := read_I1[:TAGLENGTH]
		index_i7 := read_I1[len(read_I1)-TAGLENGTH:]
		index_p5 := read_I2[:TAGLENGTH]
		index_i5 := read_I2[len(read_I2)-TAGLENGTH:]

		index := fmt.Sprintf("%s%s%s%s", index_p7, index_i7, index_p5, index_i5)
		index_1 := fmt.Sprintf("@%s:%s", index, id_I1[1:])
		index_2 := fmt.Sprintf("@%s:%s", index, id_I2[1:])

		to_write_I1 := fmt.Sprintf("%s\n%s\n%s\n%s\n", index_1, read_I1, strand_I1, qual_I1)
		to_write_I2 := fmt.Sprintf("%s\n%s\n%s\n%s\n", index_2, read_I2, strand_I2, qual_I2)
		to_write_R1 := fmt.Sprintf("%s\n%s\n%s\n%s\n", index_1, read_R1, strand_R1, qual_R1)
		to_write_R2 := fmt.Sprintf("%s\n%s\n%s\n%s\n", index_2, read_R2, strand_R2, qual_R2)

		bzip_I1.Write([]byte(to_write_I1))
		bzip_I2.Write([]byte(to_write_I2))
		bzip_R1.Write([]byte(to_write_R1))
		bzip_R2.Write([]byte(to_write_R2))

		count += 1

		if MAX_NB_READS !=0 && count > MAX_NB_READS {
			break
		}
	}

	t_end := time.Now().Second()
	fmt.Printf("demultiplexing finished in %d s\n", t_end - t_start)

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
