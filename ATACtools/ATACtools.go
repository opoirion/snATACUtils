package main

import(
	"fmt";
	"flag";
	"bufio";
	"os";
	"log";
	"time";
	"strings";
	"path";
	"sort";
	"strconv";
	// "github.com/dsnet/compress/bzip2"
	utils "ATACdemultiplex/ATACdemultiplexUtils"
	// Cbzip2 "ATACdemultiplex/cbzip2"
)




/*FILENAME ...*/
var FILENAME string
/*MAX ...*/
var MAX int
/*BZ2 ...*/
var BZ2 bool
/*GZ ...*/
var GZ bool
/*PRINTLASTLINE ...*/
var PRINTLASTLINE bool
/*PRINTLASTLINES ...*/
var PRINTLASTLINES int
/*BZ2PUREGO ...*/
var BZ2PUREGO bool
/*READTOWRITE ...*/
var READTOWRITE bool
/*COMPRESSIONMODE ...*/
var COMPRESSIONMODE int
/*GOTOLINE ...*/
var GOTOLINE int
/*SEARCHLINE ...*/
var SEARCHLINE string
/*SEP ...*/
var SEP string
/*SORTLOGS ...*/
var SORTLOGS bool
/*IGNORESORTINGCATEGORY ...*/
var IGNORESORTINGCATEGORY bool
/*IGNOREERROR ...*/
var IGNOREERROR bool


func main() {
	flag.StringVar(&FILENAME, "filename", "", "name of the file to count the lines")
	flag.IntVar(&COMPRESSIONMODE, "compressionMode", 6, "compressionMode for native bzip2 lib (1 faster -> 9 smaller) <default: 6>")
	flag.IntVar(&MAX, "max nb lines", 0, "max number of lines")
	flag.BoolVar(&BZ2, "bz2", false, "is bz2")
	flag.BoolVar(&GZ, "gz", false, "is gz")
	flag.IntVar(&GOTOLINE, "gotoline", 0, "go to line")
	flag.BoolVar(&SORTLOGS, "sortfile", false, "sort files (<key><SEP><value> file")
	flag.BoolVar(&PRINTLASTLINE, "printlastline", false, "print last line")
	flag.IntVar(&PRINTLASTLINES, "printlastlines", 0, "print last n lines")
	flag.BoolVar(&BZ2PUREGO, "bz2PureGo", false, "is bz2 using pureGo")
	flag.BoolVar(&READTOWRITE, "readtowrite", false, "read to write")
	flag.BoolVar(&IGNOREERROR, "ignoreerror", false, "ignore error and continue")
	flag.BoolVar(&IGNORESORTINGCATEGORY, "ignore_sorting_category", false, "ignore file cateogry (identified by #) when sorting")
	flag.StringVar(&SEARCHLINE, "search_in_line", "", "search specific motifs in line")
	flag.StringVar(&SEP, "delimiter", "\t", "delimiter used to split and sort the log file (default \t)")
	flag.Parse()
	fmt.Printf("%s\n", FILENAME)
	tStart := time.Now()

	var nbLines int

	switch  {
	case SORTLOGS:
		sortLogfile(FILENAME, SEP)
		nbLines = 0
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

	tDiff := time.Now().Sub(tStart)

	if nbLines > 0 {
		fmt.Printf("number of lines: %d\n", nbLines)
	}
	fmt.Printf("time: %f s\n", tDiff.Seconds())

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
	writer := utils.ReturnWriterForBzipfile("tmp.bz2", COMPRESSIONMODE)

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

//  ...
func sortLogfile(filename string, separator string)  {
	file, err := os.Open(filename)
	defer file.Close()
	check(err)
	scanner := bufio.NewScanner(file)

	ext := path.Ext(filename)
	outfname := fmt.Sprintf("%s_sorted%s", strings.TrimSuffix(filename, ext), ext)
	outfile, err := os.Create(outfname)
	defer outfile.Close()

	defer os.Rename(outfname, filename)

	check(err)
	pl := utils.PairList{}

	lineNb := -1
	buff := 0

	for scanner.Scan() {
		line := scanner.Text()
		lineNb += 1

		if len(line) == 0 || line[0] == '#' || line[0] == '\n'  {
			if IGNORESORTINGCATEGORY{
				continue
			}

			outfile.WriteString(fmt.Sprintf("%s\n", line))
			if len(pl) == 0  {
				continue
			}
			sort.Sort(sort.Reverse(pl))

			for _, el := range pl {
				outfile.WriteString(fmt.Sprintf("%s%s%d\n", el.Key, SEP, el.Value))

				buff += 1
				if buff > 100000{
					outfile.Sync()
					buff = 0
				}
			}

			pl = utils.PairList{}
			continue
		}
		split := strings.Split(line, SEP)
		valueField := split[len(split)-1]
		key := strings.Join(split[:len(split)-1], SEP)
		value, err := strconv.Atoi(valueField)

		if err != nil {
			fmt.Printf("value field %s from: %s at line %d not conform!\n",
				valueField, line, lineNb)

			if !IGNOREERROR {
				check(err)
			}
		}

		pl = append(pl, utils.Pair{Key:key, Value:value})
	}

	if len(pl) == 0 {
		return
	}

	sort.Sort(sort.Reverse(pl))

	for _, el := range pl {
		outfile.WriteString(fmt.Sprintf("%s%s%d\n", el.Key, SEP, el.Value))
		buff += 1
		if buff > 100000{
			outfile.Sync()
			buff = 0
		}
	}
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
