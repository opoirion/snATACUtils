package main

import(
	"fmt";
	"flag";
	"bufio";
	"os";
	"log";
	"time";
	"strings";
	// "github.com/dsnet/compress/bzip2"
	utils "ATACdemultiplex/ATACdemultiplexUtils"
	// Cbzip2 "ATACdemultiplex/cbzip2"
)




var FILENAME string
var MAX int
var BZ2 bool
var GZ bool
var PRINTLASTLINE bool
var PRINTLASTLINES int
var BZ2PUREGO bool
var READTOWRITE bool
var COMPRESSION_MODE int
var GOTOLINE int
var SEARCHLINE string

func main() {
	flag.StringVar(&FILENAME, "filename", "", "name of the file to count the lines")
	flag.IntVar(&COMPRESSION_MODE, "compressionMode", 6, "compressionMode for native bzip2 lib (1 faster -> 9 smaller) <default: 6>")
	flag.IntVar(&MAX, "max nb lines", 0, "max number of lines")
	flag.BoolVar(&BZ2, "bz2", false, "is bz2")
	flag.BoolVar(&GZ, "gz", false, "is gz")
	flag.IntVar(&GOTOLINE, "gotoline", 0, "go to line")
	flag.BoolVar(&PRINTLASTLINE, "printlastline", false, "print last line")
	flag.IntVar(&PRINTLASTLINES, "printlastlines", 0, "print last n lines")
	flag.BoolVar(&BZ2PUREGO, "bz2PureGo", false, "is bz2 using pureGo")
	flag.BoolVar(&READTOWRITE, "readtowrite", false, "read to write")
	flag.StringVar(&SEARCHLINE, "search_in_line", "", "search specific motifs in line")
	flag.Parse()
	fmt.Printf("%s\n", FILENAME)
	t_start := time.Now()

	var nbLines int

	switch  {
	case BZ2:
		nbLines = countLineBz2(FILENAME)
	case GZ:
		nbLines = countLineGz(FILENAME)
	case BZ2PUREGO:
		nbLines = countLineBz2PureGo(FILENAME, 0)
	case READTOWRITE:
		nbLines = readtowrite(FILENAME)
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
	return processScanner(scanner)
}

// readtowrite ...
func readtowrite(filename string) int  {
	scanner, file := utils.ReturnReaderForBzipfile(filename, 0)
	defer file.Close()
	writer := utils.ReturnWriterForBzipfile("tmp.bz2", COMPRESSION_MODE)

	nbLines := 0

	for scanner.Scan() {
		nbLines++
		line := scanner.Text()
		if len(line) > 6 {

			line = fmt.Sprintf(">%s%s%s", line[:4], line[:5], strings.Split(line,"@" )[0])

		}

		writer.Write([]byte(line))
	}

	utils.ExceCmd("rm tmp.bz2")

	return nbLines

}

func countLineBz2(filename string) int {

	scanner, file := utils.ReturnReaderForBzipfile(filename, 0)
	defer file.Close()

	return processScanner(scanner)
}

func processScanner(scanner * bufio.Scanner) (nbLines int) {
	nbLines = 0
	for scanner.Scan() {
		nbLines++

		if (GOTOLINE > 0) && (PRINTLASTLINES > 0) && (GOTOLINE - nbLines < PRINTLASTLINES)  {
			line := scanner.Text()
			fmt.Printf("last lines: %s\n", line)
		}

		if GOTOLINE > 0 && nbLines > GOTOLINE {
			break
		}

		if SEARCHLINE != "" {
			line := scanner.Text()
			if strings.Contains(line, SEARCHLINE) {
				fmt.Printf("line nb: %d\t %s found in %s\n", nbLines, SEARCHLINE, line)
			}
		}
	}

	if PRINTLASTLINE {
		line := scanner.Text()
		fmt.Printf("last line: %s\n", line)
	}

	return nbLines

}

func countLineGz(filename string) int {

	scanner, file := utils.ReturnReaderForGzipfile(filename, 0)
	defer file.Close()

	return processScanner(scanner)
}

func countLineBz2PureGo(filename string, pos int) int {

	scanner, file := utils.ReturnReaderForBzipfilePureGo(filename, pos)
	defer file.Close()

	return processScanner(scanner)
}


func check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
