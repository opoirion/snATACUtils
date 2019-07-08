package main

import (
	"os"
	"bufio"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"path"
	"fmt"
)


func countLine(fname string) int {
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

	for scanner.Scan() {
		i++
		if i>nbLines {break}
	}

	seek, err := fileScanner.Seek(0, 1)
	utils.Check(err)
	fmt.Printf("estimated number of lines: %d\n", (int(size) / int(seek)) * i)

	return (int(size) / int(seek)) * i

}
