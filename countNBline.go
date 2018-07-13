package main

import (
	"os"
	"log"
	"bufio"
	originalbzip2  "compress/bzip2"
	utils "ATACdemultiplex/ATACdemultiplexUtils"
	"strings"
	"path"
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

	nbLines := 50000
	i := 0

	var scanner * bufio.Scanner
	var fileScanner * os.File

	switch ext {
	case ".bz2":
		scanner, fileScanner = utils.ReturnReaderForBzipfile(fname, 0)
	case ".gz":
		scanner, fileScanner = utils.ReturnReaderForGzipfile(fname, 0)
	default:
		panic(fmt.Sprintf("%s not a valid extension for %s!", ext, fname))
	}


	defer fileScanner.Close()
	// defer writer.Close()

	for scanner.Scan() {
		i++
		if i>nbLines {break}
	}

	seek, err := fileScanner.Seek(0, 1)
	fmt.Printf("estimated number of lines: %d\n", (int(size) / int(seek)) * i)

	return (int(size) / int(seek)) * i

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
