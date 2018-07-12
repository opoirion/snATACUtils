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


/*FILENAME ...*/
var FILENAME string
/*MAX ...*/
var MAX int

func countLine(fname string, compressionMode int) int {
	fStat, err := os.Stat(fname)
	Check(err)

	ext := path.Ext(fname)

	size := fStat.Size()

	nbLines := 4000

	var scanner * bufio.Scanner
	var fileScanner * os.File
	var writer io.WriteCloser
	var factor int

	switch ext {
	case ".bz2":
		scanner, fileScanner = utils.ReturnReaderForBzipfile(fname, 0)
		writer = utils.ReturnWriterForBzipfile(fmt.Sprintf("tmp%s", ext), compressionMode)
		factor = 1
	case ".gz":
		scanner, fileScanner = utils.ReturnReaderForGzipfile(fname, 0)
		writer = utils.ReturnWriterForGzipFile(fmt.Sprintf("tmp%s", ext))
		factor = 2
	default:
		panic(fmt.Sprintf("%s not a valid extension for %s!", ext, fname))
	}


	defer fileScanner.Close()
	defer writer.Close()

	for i:=0;i<nbLines;i++ {
		scanner.Scan()
		writer.Write([]byte(scanner.Text()))
	}
	writer.Close()

	fStatWriter, err := os.Stat(fmt.Sprintf("tmp%s", ext))
	Check(err)

	sizeWriter := fStatWriter.Size()
	utils.ExceCmd(fmt.Sprintf("rm tmp%s", ext))

	return (int(size) / int(sizeWriter)) * int(nbLines) / factor

}

func guessPosToGo(fname string, lineToGo int) int {
	fileOpen, err := os.OpenFile(fname, 0, 0)

	if err != nil {
		log.Fatal(err)
	}

	readerOs := bufio.NewReader(fileOpen)
	readerBzip := originalbzip2.NewReader(readerOs)

	chunk := make([]byte, 10000)
	readerBzip.Read(chunk)
	nbLines := strings.Count(string(chunk), "\n")

	return lineToGo * 10000 / nbLines

}
