package main

import (
	"os"
	"fmt"
	"log"
	"bufio"
	originalbzip2  "compress/bzip2"
	utils "ATACdemultiplex/ATACdemultiplexUtils"
	"strings"
)


var FILENAME string
var MAX int

func countLine(fname string, compressionMode int) int {
	f_stat, err := os.Stat(fname)
	Check(err)

	size := f_stat.Size()
	scanner, file_scanner := utils.ReturnReaderForBzipfile(fname, 0)
	defer file_scanner.Close()
	writer := utils.ReturnWriterForBzipfile("tmp.bz2", compressionMode)
	defer writer.Close()

	for i:=0;i<2000;i++ {
		scanner.Scan()
		writer.Write([]byte(scanner.Text()))
	}
	writer.Close()

	f_stat_writer, err := os.Stat("tmp.bz2")
	Check(err)

	size_writer := f_stat_writer.Size()
	utils.ExceCmd("rm tmp.bz2")

	return (int(size) / int(size_writer)) * int(2000)

}

func guessPosToGo(fname string, lineToGo int) int {
	file_open, err := os.OpenFile(fname, 0, 0)

	if err != nil {
		log.Fatal(err)
	}

	reader_os := bufio.NewReader(file_open)
	reader_bzip := originalbzip2.NewReader(reader_os)

	chunk := make([]byte, 10000)
	reader_bzip.Read(chunk)
	nbLines := strings.Count(string(chunk), "\n")

	fmt.Printf("number of lines for %s: %d\n", fname, nbLines)

	return lineToGo * 10000 / nbLines

}
