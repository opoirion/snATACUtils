package main

import(
	"fmt";
	"flag";
	"bufio";
	"os";
	"log";
	"time";
	"github.com/dsnet/compress/bzip2"
)




var FILENAME string
var MAX int
var BZ2 bool

func main() {
	flag.StringVar(&FILENAME, "filename", "", "name of the file to count the lines")
	flag.IntVar(&MAX, "max nb lines", 0, "max number of lines")
	flag.BoolVar(&BZ2, "bz2", false, "is bz2")
	flag.Parse()
	fmt.Printf("%s\n", FILENAME)
	t_start := time.Now()

	var nbLines int

	switch  {
	case BZ2:
		nbLines = countLineBz2(FILENAME)
	default:
		nbLines = countLine(FILENAME)

	}

	t_diff := time.Now().Sub(t_start)

	fmt.Printf("number of lines: %d\n", nbLines)
	fmt.Printf("time: %f s\n", t_diff.Seconds())

}


func countLine(filename string) int {
	file, err := os.Open(filename)
	check(err)

	scanner := bufio.NewScanner(file)
	nbLines := 0

	for scanner.Scan() {
		nbLines++
	}

	return nbLines
}

func countLineBz2(filename string) int {

	scanner := return_reader_for_bzipfile(filename)
	nbLines := 0

	for scanner.Scan() {
		nbLines++
	}

	return nbLines
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


func check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
