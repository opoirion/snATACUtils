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
	"sync";
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
/*CREATEREFFASTQ ...*/
var CREATEREFFASTQ bool
/*REFBARCODELIST ...*/
var REFBARCODELIST string
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
/*COMPLSTRATEGY ...*/
var COMPLSTRATEGY string
/*CREATEBARCODEDICT ...*/
var CREATEBARCODEDICT bool
/*FILENAMES ...*/
var FILENAMES utils.ArrayFlags
/*TAG ...*/
var TAG string


func main() {
	flag.StringVar(&FILENAME, "filename", "", "name of the file(s) multiple files should be into \"")
	flag.Var(&FILENAMES, "filenames", "name of the files to use")
	flag.IntVar(&COMPRESSIONMODE, "compressionMode", 6, "compressionMode for native bzip2 lib (1 faster -> 9 smaller) <default: 6>")
	flag.IntVar(&MAX, "max_nb_lines", 0, "max number of lines")
	flag.BoolVar(&BZ2, "bz2", false, "is bz2")
	flag.BoolVar(&GZ, "gz", false, "is gz")
	flag.BoolVar(&CREATEREFFASTQ, "create_ref_fastq", false, "create a ref FASTQ file using a reference barcode list")
	flag.StringVar(&REFBARCODELIST, "ref_barcode_list", "", "file containing the reference barcodes (one per line)")
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
	flag.StringVar(&TAG, "tag", "", "tag used when creating a reference fastq file to tag all the reads (default \"\")")
	flag.StringVar(&COMPLSTRATEGY, "compl_strategy", "split_10_compl_second", "Strategy to use when writing the complement of a fastq file (default split_10_compl_second)")
	flag.BoolVar(&CREATEBARCODEDICT, "create_barcode_dict", false, "create a barcode key / value count file")
	flag.Parse()
	fmt.Printf("input file(s): %s\n", FILENAME)
	tStart := time.Now()

	var nbLines int

	switch  {
	case WRITECOMPL:
		writeComplement(FILENAME, COMPLSTRATEGY)
	case CREATEREFFASTQ && len(FILENAMES) > 0:
		var waiting sync.WaitGroup
		waiting.Add(len(FILENAMES))

		for _, filename := range(FILENAMES){
			go extractFASTQreadsPerBarcodes(filename, REFBARCODELIST, &waiting)
		}
		waiting.Wait()

	case CREATEREFFASTQ:
		var waiting sync.WaitGroup
		extractFASTQreadsPerBarcodes(FILENAME, REFBARCODELIST, &waiting)
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

/*extractFASTQreadsPerBarcodes create a new FASTQ file using a barcode name */
func extractFASTQreadsPerBarcodes(filename string, barcodefilename string, waiting *sync.WaitGroup) {

	var barcode string
	var refbarcodes = make(map[string]bool)
	var readHeader = make([]string, 2, 2)
	defer waiting.Done()

	ext := path.Ext(filename)
	ext2 := path.Ext(filename[:len(filename) - len(ext)])

	outfile := fmt.Sprintf("%s.reference%s%s",
		filename[:len(filename) - len(ext) -len(ext2)], ext2, ext)

	scanner, file := utils.ReturnReader(filename, 0, false)
	defer file.Close()
	writer := utils.ReturnWriter(outfile, COMPRESSIONMODE, false)
	defer writer.Close()

	barcodefile, err := os.Open(barcodefilename)
	check(err)
	defer barcodefile.Close()
	bscanner := bufio.NewScanner(barcodefile)

	for bscanner.Scan() {
		line := bscanner.Text()
		refbarcodes[line] = true
	}

	isfour := 0
	nbLine := 0
	tocopy := false

	for scanner.Scan() {
		line := scanner.Text()
		nbLine++

		if (isfour == 0) {
			readHeader = strings.SplitN(line, ":", 2)
			barcode = readHeader[0][1:]
			if refbarcodes[barcode] {
				tocopy = true
				nbLine++
			}
		}

		if tocopy {
			if (isfour == 0) && (len(TAG) > 0) {
				writer.Write([]byte(fmt.Sprintf("@%s%s:%s\n", barcode, TAG, readHeader[1])))
			} else {
				writer.Write([]byte(fmt.Sprintf("%s\n", line)))
			}
		}

		isfour++

		if isfour == 4 {
			tocopy = false
			isfour = 0
		}
	}
	fmt.Printf("nb barcodes extracted: %d\n", nbLine)
	fmt.Printf("file created: %s\n", outfile)
}


func writeComplement(filename string, complStrategy string) (nbLines int) {

	ext := path.Ext(filename)
	ext2 := path.Ext(filename[:len(filename) - len(ext)])

	outfile := fmt.Sprintf("%s.compl%s%s",
		filename[:len(filename) - len(ext) -len(ext2)], ext2, ext)

	scanner, file := utils.ReturnReader(filename, 0, false)
	defer file.Close()
	writer := utils.ReturnWriter(outfile, COMPRESSIONMODE, false)
	defer writer.Close()

	cycle := 0

	var newID string

	for scanner.Scan() {
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

			switch COMPLSTRATEGY {
			case "split_10_compl_second":
				split1 := barcode[:10]
				split2 := barcode[10:]

				newID = fmt.Sprintf("@%s%s:%s\n",
					split1, returnComp(split2), rest)

			case "split_10_compl_first":
				split1 := barcode[:10]
				split2 := barcode[10:]

				newID = fmt.Sprintf("@%s%s:%s\n",
					returnComp(split1), split2, rest)
			default:
				panic(fmt.Sprintf("wrong strategy used!: %s\n", COMPLSTRATEGY))
			}

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

	if CREATEBARCODEDICT {
		defer createIndexCountFile(outfile)
	}

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
		defer sortLogfile(outfilename, SEP)
	}

	fmt.Printf("output file: %s\n", outfilename)

	return nbLines
}

func returnComp(s string) (comp []byte) {
	length := len(s)
	comp = make([]byte, length)

	for pos, char := range s {
		if _, ok := COMPDICT[byte(char)];!ok {
			panic(
				fmt.Sprintf(
				"error with substring when finding the complementary: %s elem not in compdict\n", s))
		}

		comp[length - pos - 1] = COMPDICT[byte(char)]
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
		defer sortLogfile(outfname, SEP)
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
