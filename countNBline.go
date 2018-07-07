package main

import (
	"os"
	"log"
	"bufio"
	originalbzip2  "compress/bzip2"
	utils "ATACdemultiplex/ATACdemultiplexUtils"
	"strings"
	"path"
	"io"
	"fmt"
)


var FILENAME string
var MAX int

func countLine(fname string, compressionMode int) int {
	f_stat, err := os.Stat(fname)
	Check(err)

	ext := path.Ext(fname)

	size := f_stat.Size()

	nbLines := 4000

	var scanner * bufio.Scanner
	var file_scanner * os.File
	var writer io.WriteCloser
	var factor int

	switch ext {
	case ".bz2":
		scanner, file_scanner = utils.ReturnReaderForBzipfile(fname, 0)
		writer = utils.ReturnWriterForBzipfile(fmt.Sprintf("tmp%s", ext), compressionMode)
		factor = 1
	case ".gz":
		scanner, file_scanner = utils.ReturnReaderForGzipfile(fname, 0)
		writer = utils.ReturnWriterForGzipFile(fmt.Sprintf("tmp%s", ext))
		factor = 1
	default:
		panic(fmt.Sprintf("%s not a valid extension for %s!", ext, fname))
	}


	defer file_scanner.Close()
	defer writer.Close()

	for i:=0;i<nbLines;i++ {
		scanner.Scan()
		writer.Write([]byte(scanner.Text()))
	}
	writer.Close()

	f_stat_writer, err := os.Stat(fmt.Sprintf("tmp%s", ext))
	Check(err)

	size_writer := f_stat_writer.Size()
	utils.ExceCmd(fmt.Sprintf("rm tmp%s", ext))

	return (int(size) / int(size_writer)) * int(nbLines) / factor

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

	return lineToGo * 10000 / nbLines

}
