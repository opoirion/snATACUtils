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
/*OUTFILE ...*/
var OUTFILE string
/*MERGE ...*/
var MERGE bool
/*WRITECOMPL ...*/
var WRITECOMPL bool
/*COMPDICT ...*/
var COMPDICT = map[byte]byte {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
/*SPLITBARCODEAFTER ...*/
var SPLITBARCODEAFTER int
/*CREATEBARCODEDICT ...*/
var CREATEBARCODEDICT bool


func main() {
	flag.StringVar(&FILENAME, "filename", "", "name of the file(s) multiple files should be into \"")
	flag.IntVar(&COMPRESSIONMODE, "compressionMode", 6, "compressionMode for native bzip2 lib (1 faster -> 9 smaller) <default: 6>")
	flag.IntVar(&MAX, "max_nb_lines", 0, "max number of lines")
	flag.BoolVar(&BZ2, "bz2", false, "is bz2")
	flag.BoolVar(&GZ, "gz", false, "is gz")
	flag.BoolVar(&MERGE, "merge", false, "merge input log files together")
	flag.IntVar(&GOTOLINE, "gotoline", 0, "go to line")
	flag.BoolVar(&SORTLOGS, "sortfile", false, "sort files (<key><SEP><value> file")
	flag.BoolVar(&PRINTLASTLINE, "printlastline", false, "print last line")
	flag.IntVar(&PRINTLASTLINES, "printlastlines", 0, "print last n lines")
	flag.BoolVar(&BZ2PUREGO, "bz2PureGo", false, "is bz2 using pureGo")
	flag.BoolVar(&READTOWRITE, "readtowrite", false, "read to write")
	flag.BoolVar(&IGNOREERROR, "ignoreerror", false, "ignore error and continue")
	flag.BoolVar(&IGNORESORTINGCATEGORY, "ignore_sorting_category", false, "ignore file cateogry (identified by #) when sorting")
	flag.StringVar(&SEARCHLINE, "search_in_line", "", "search specific motifs in line")
	flag.StringVar(&OUTFILE, "output", "", "file name of the output")
	flag.BoolVar(&WRITECOMPL, "write_compl", false, "write the barcode complement of a fastq files")
	flag.StringVar(&SEP, "delimiter", "\t", "delimiter used to split and sort the log file (default \t)")
	flag.IntVar(&SPLITBARCODEAFTER, "splitbarcodeafter", 10, "when writing the complement of a fastq file, split the barcode ID after N nucleotides (default 10)")
	flag.BoolVar(&CREATEBARCODEDICT, "create_barcode_dict", false, "create a barcode key / value count file")
	flag.Parse()
	fmt.Printf("input file(s): %s\n", FILENAME)
	tStart := time.Now()

	var nbLines int

	switch  {
	case WRITECOMPL:
		writeComplement(FILENAME, SPLITBARCODEAFTER)
	case CREATEBARCODEDICT:
		createIndexCountFile(FILENAME)
	case MERGE:
		mergeLogFiles(strings.Split(FILENAME, " "), OUTFILE)
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


func writeComplement(filename string, splitBarcodeAfter int) (nbLines int) {

	ext := path.Ext(filename)
	ext2 := path.Ext(filename[:len(filename) - len(ext)])

	outfile := fmt.Sprintf("%s.compl%s%s",
		filename[:len(filename) - len(ext) -len(ext2)], ext2, ext)

	scanner, file := utils.ReturnReader(filename, 0, false)
	defer file.Close()
	writer := utils.ReturnWriter(outfile, COMPRESSIONMODE, false)
	defer writer.Close()

	cycle := 0

	for scanner .Scan() {
		nbLines++
		line := scanner.Text()

		switch cycle {
		case 0:

			split := strings.SplitN(line, ":", 2)

			barcode := split[0][1:]
			rest := split[1]

			if split[0][0] != '@' {
				panic(fmt.Sprintf("idLine not conform:%s at pos:%d\n", line, nbLines))
			}

			split1 := barcode[:splitBarcodeAfter]
			split2 := barcode[splitBarcodeAfter:]

			newID := fmt.Sprintf("@%s%s:%s\n",
				returnComp(split1), returnComp(split2), rest)

			writer.Write([]byte(newID))
		default:
			writer.Write([]byte(line + "\n"))

		}

		cycle++

		if cycle == 4 {
			cycle = 0
		}

		if (MAX > 0) && (MAX < nbLines) {
			break
		}

	}

	fmt.Printf("output file: %s\n", outfile)

	return nbLines
}

func createIndexCountFile(filename string) (nbLines int) {

	ext := path.Ext(filename)
	ext2 := path.Ext(filename[:len(filename) - len(ext)])

	outfilename := fmt.Sprintf("%s.barcodeCounts",
		filename[:len(filename) - len(ext) -len(ext2)])

	scanner, file := utils.ReturnReader(filename, 0, false)
	defer file.Close()
	outfile, err := os.Create(outfilename)
	check(err)

	defer outfile.Close()
	dict := map[string]int {}

	cycle := 0

	for scanner .Scan() {
		nbLines++
		line := scanner.Text()

		switch cycle {
		case 0:

			split := strings.SplitN(line, ":", 2)

			barcode := split[0][1:]

			dict[barcode]++
		}

		cycle++

		if cycle == 4 {
			cycle = 0
		}

		if (MAX > 0) && (MAX < nbLines) {
			break
		}
	}

	for key, value := range dict {
		outfile.WriteString(fmt.Sprintf("%s%s%d\n", key, SEP, value))
	}

	if SORTLOGS {
		outfile.Close()
		sortLogfile(outfilename, SEP)
	}

	fmt.Printf("output file: %s\n", outfilename)

	return nbLines
}

func returnComp(s string) (comp []byte) {
	comp = []byte(s)

	for pos, char := range comp {
		if _, ok := COMPDICT[char];!ok {
			panic(
				fmt.Sprintf(
				"error with substring when finding the complementary: %s elem not in compdict\n", s))
		}

		comp[pos] = COMPDICT[char]
	}

	return comp
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

		key, value := splitLine(line, lineNb)
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


func splitLine(line string, lineNb int) (key string, value int){
	split := strings.Split(line, SEP)
	valueField := split[len(split)-1]
	key = strings.Join(split[:len(split)-1], SEP)
	value, err := strconv.Atoi(valueField)

	if err != nil {
		fmt.Printf("value field %s from: %s at line nb %d not conform!\n",
			valueField, line, lineNb)

		if !IGNOREERROR {
			check(err)
		}
	}

	return key, value
}

func mergeLogFiles(filenames []string, outfname string) {
	fmt.Printf("output file: %s\n", OUTFILE)
	if outfname == "" {
		panic("option -output should be non null!")
	}

	outfile, err := os.Create(outfname)
	check(err)
	defer outfile.Close()

	var key string
	var value int

	dict := map[string]int{}

	for _, filename := range filenames {

		if filename == "" {
			continue
		}

		file, err := os.Open(filename)
		if err != nil {
			fmt.Printf("filename: %s presents error: %s continue...\n", filename, err)
			continue
		}

		defer file.Close()

		scanner := bufio.NewScanner(file)

		for scanner.Scan() {
			line := scanner.Text()

			if len(line) == 0 || line[0] == '#' || line[0] == '\n'  {
				continue
			}
			key, value = splitLine(line, 0)
			dict[key] += value
		}
	}

	buff := 0

	for key, value := range dict {
		outfile.WriteString(fmt.Sprintf("%s%s%d\n", key, SEP, value))

		buff++
		if buff > 100000{
			outfile.Sync()
			buff = 0
		}
	}

	if SORTLOGS {
		outfile.Close()
		sortLogfile(outfname, SEP)
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
