/* Suite of functions dedicated to pre/post process files related to snATAC pipeline */

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
	"sync";
	"strconv";
	"bytes"
	utils "ATACdemultiplex/ATACdemultiplexUtils"
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
/*CREATEREFBEDFILE ...*/
var CREATEREFBEDFILE bool
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
/*OUTTAG ...*/
var OUTTAG string
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
/*CICEROPROCESSING ...*/
var CICEROPROCESSING bool


func main() {
	flag.StringVar(&FILENAME, "filename", "", "name of the file(s) multiple files should be into \"")
	flag.Var(&FILENAMES, "filenames", "name of the files to use")
	flag.IntVar(&COMPRESSIONMODE, "compressionMode", 6, "compressionMode for native bzip2 lib (1 faster -> 9 smaller) <default: 6>")
	flag.IntVar(&MAX, "max_nb_lines", 0, "max number of lines")
	flag.BoolVar(&BZ2, "bz2", false, "is bz2")
	flag.BoolVar(&GZ, "gz", false, "is gz")
	flag.BoolVar(&CICEROPROCESSING, "bed_to_cicero", false,
		`format a bed to cicero input (ex: chr1\t1215523\t1216200\tcellID -> chr1_1215523_121620,tcellID,1)
            USAGE: ATACtools -bed_to_cicero -filename <bedfile> (-filenames <bedfile2> -filenames <bedfile3> ...  -ignoreerror)`)
	flag.BoolVar(&CREATEREFFASTQ, "create_ref_fastq", false, `create a ref FASTQ file using a reference barcode list
            USAGE: ATACtools -create_ref_bed -filename <fname> (-ref_barcode_list <fname> -tag <string>)`)
	flag.BoolVar(&CREATEREFBEDFILE, "create_ref_bed", false, `create a ref bed file using a reference barcode list
            USAGE: ATACtools -create_ref_fastq -filnames <fname1> -filnames <fname2> ... (-ref_barcode_list <fname> -tag <string>)`)
	flag.StringVar(&REFBARCODELIST, "ref_barcode_list", "", "file containing the reference barcodes (one per line)")
	flag.BoolVar(&MERGE, "merge", false, `merge input log files together
            USAGE: ATACtools -merge -filenames <fname1> -filenames <fname2> -filenames <fname3> ...  (-sortfile -delimiter "<string>" -ignoreerror -ignore_sorting_category)`)
	flag.IntVar(&GOTOLINE, "gotoline", 0, "go to line")
	flag.BoolVar(&SORTLOGS, "sortfile", false,
		`sort files (<key><SEP><value> file
            USAGE: ATACtools -sortfile -filename <fname> (-delimiter <string> -ignoreerror -ignore_sorting_category)`)
	flag.BoolVar(&PRINTLASTLINE, "printlastline", false, "print last line")
	flag.IntVar(&PRINTLASTLINES, "printlastlines", 0, "print last n lines")
	flag.BoolVar(&BZ2PUREGO, "bz2PureGo", false, "is bz2 using pureGo")
	flag.BoolVar(&READTOWRITE, "readtowrite", false, "read to write")
	flag.BoolVar(&IGNOREERROR, "ignoreerror", false, "ignore error and continue")
	flag.BoolVar(&IGNORESORTINGCATEGORY, "ignore_sorting_category", false, "ignore file cateogry (identified by #) when sorting")
	flag.StringVar(&SEARCHLINE, "search_in_line", "", "search specific motifs in line")
	flag.StringVar(&OUTFILE, "output", "", "file name of the output")
	flag.StringVar(&OUTTAG, "output_tag", "reference", "particule to annotate the output file name")
	flag.BoolVar(&WRITECOMPL, "write_compl", false, `write the barcode complement of a fastq files
            USAGE: ATACtools -write_compl <fastq_file> (-compl_strategy <"split_10_compl_second"/"split_10_compl_first"> -tag <string>)`)
	flag.StringVar(&SEP, "delimiter", "\t", "delimiter used to split and sort the log file (default \t)")
	flag.StringVar(&TAG, "tag", "", "tag used when creating a reference fastq file to tag all the reads (default \"\")")
	flag.StringVar(&COMPLSTRATEGY, "compl_strategy", "split_10_compl_second", `Strategy to use when writing the complement of a fastq file (default split_10_compl_second: split after 10 bases and complementary only second)`)
	flag.BoolVar(&CREATEBARCODEDICT, "create_barcode_dict", false, `create a barcode key / value count file
            USAGE: ATACtools -create_barcode_list -filename <fname> (-sortfile -delimiter <string>)`)
	flag.Parse()


	if FILENAME != "" {
		fmt.Printf("input file(s): %s\n", FILENAME)
	}

	for _, filename := range FILENAMES {
		fmt.Printf("input file(s): %s\n", filename)
	}

	tStart := time.Now()

	var nbLines int

	switch  {
	case CICEROPROCESSING:
		bedFilestoCiceroInput(FILENAMES, OUTFILE)

	case WRITECOMPL:
		writeComplement(FILENAME, COMPLSTRATEGY)

	case CREATEREFBEDFILE:
		extractBEDreadsPerBarcodes(FILENAME, REFBARCODELIST)

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
		mergeLogFiles(FILENAMES, OUTFILE)
	case SORTLOGS:
		utils.SortLogfile(FILENAME, SEP, "", IGNORESORTINGCATEGORY, IGNOREERROR)
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
/*bedFilestoCiceroInput ...*/
func bedFilestoCiceroInput(filenames []string, outfile string) {
	var ext string
	var finalOutfile string

	if FILENAME != "" {
		filenames = append(filenames, FILENAME)
	}

	if len(filenames) == 0 {
		log.Fatal("at least one input file (option -filename(s) <string>) must be filled!")
	}

	if outfile == "" {
		ext = path.Ext(filenames[0])
		finalOutfile = fmt.Sprintf("%s.cicero_input%s", filenames[0][:len(filenames[0])-len(ext)], ext)
	} else {

		finalOutfile = outfile
	}

	var waiting sync.WaitGroup
	waiting.Add(len(filenames))

	outfiles := []string{}

	tStart := time.Now()
	ext = path.Ext(finalOutfile)

	for i, file := range filenames {
		outfile := fmt.Sprintf("%s.index_%d%s", finalOutfile, i, ext)
		go bedFiletoCiceroInputOnethread(file, outfile, REFBARCODELIST, &waiting, i==0)
		outfiles = append(outfiles, outfile)
	}

	waiting.Wait()

	cmd := fmt.Sprintf("cat %s > %s", strings.Join(outfiles, " "), finalOutfile)
	utils.ExceCmd(cmd)
	cmd = fmt.Sprintf("rm %s", strings.Join(outfiles, " "))
	utils.ExceCmd(cmd)

	tDiff := time.Now().Sub(tStart)
	fmt.Printf("file:%s created in %f s\n", finalOutfile, tDiff.Seconds())
}

/*bedFiletoCiceroInputOnethread ...*/
func bedFiletoCiceroInputOnethread(filename string, outfile string, barcodefilename string,
	waiting * sync.WaitGroup, writeHeader bool) {
	defer waiting.Done()
	var buffer bytes.Buffer
	var split = make([]string, 4, 4)
	var nbLine int
	var barcodeIndex map[string]bool
	var checkBarcode bool

	if barcodefilename != "" {
		checkBarcode = true
		barcodeIndex = loadCellIDDict(barcodefilename)
	}

	scanner, file := utils.ReturnReader(filename, 0, false)
	defer file.Close()
	writer := utils.ReturnWriter(outfile, COMPRESSIONMODE, false)
	defer writer.Close()

	if writeHeader {
		writer.Write([]byte("Peak,cell,Count\n"))
	}

	for scanner.Scan() {
		line := scanner.Text()
		nbLine++

		split = strings.Split(line, "\t")

		if checkBarcode && !barcodeIndex[split[3]] {
			continue
		}

		buffer.WriteString(split[0])
		buffer.WriteRune('_')
		buffer.WriteString(split[1])
		buffer.WriteRune('_')
		buffer.WriteString(split[2])
		buffer.WriteRune(',')
		buffer.WriteString(split[3])
		buffer.WriteRune(',')
		buffer.WriteRune('1')
		buffer.WriteRune('\n')
		writer.Write(buffer.Bytes())
		buffer.Reset()
	}
}

/*extractBEDreadsPerBarcodes create a new FASTQ file using a barcode name */
func extractBEDreadsPerBarcodes(filename string, barcodefilename string) {

	var barcode string
	var refbarcodes map[string]bool
	var split = make([]string, 4, 4)
	var notCheckBarcode = true
	var buffer bytes.Buffer

	ext := path.Ext(filename)
	ext2 := path.Ext(filename[:len(filename) - len(ext)])

	outfile := fmt.Sprintf("%s.%s%s%s",
		filename[:len(filename) - len(ext) -len(ext2)], OUTTAG, ext2, ext)

	scanner, file := utils.ReturnReader(filename, 0, false)
	defer file.Close()
	writer := utils.ReturnWriter(outfile, COMPRESSIONMODE, false)
	defer writer.Close()

	if barcodefilename != "" {
		notCheckBarcode = false
		refbarcodes = loadCellIDDict(barcodefilename)
	}

	nbLine := 0

	for scanner.Scan() {
		line := scanner.Text()
		nbLine++

		split = strings.Split(line, "\t")
		barcode = split[3]


		if refbarcodes[barcode] || notCheckBarcode {

			if (len(TAG) > 0) {
				buffer.WriteString(split[0])
				buffer.WriteRune('\t')
				buffer.WriteString(split[1])
				buffer.WriteRune('\t')
				buffer.WriteString(split[2])
				buffer.WriteRune('\t')
				buffer.WriteString(barcode)
				buffer.WriteString(TAG)
				buffer.WriteRune('\n')

			} else {
				buffer.WriteString(line)
				buffer.WriteRune('\n')
			}

			writer.Write(buffer.Bytes())
			buffer.Reset()

		}
	}
	fmt.Printf("nb barcodes extracted: %d\n", nbLine)
	fmt.Printf("file created: %s\n", outfile)
}

/*extractFASTQreadsPerBarcodes create a new FASTQ file using a barcode name */
func extractFASTQreadsPerBarcodes(filename string, barcodefilename string, waiting *sync.WaitGroup) {

	var barcode string
	var refbarcodes = make(map[string]bool)
	var readHeader = make([]string, 2, 2)
	defer waiting.Done()
	var buffer bytes.Buffer

	ext := path.Ext(filename)
	ext2 := path.Ext(filename[:len(filename) - len(ext)])

	outfile := fmt.Sprintf("%s.%s%s%s",
		filename[:len(filename) - len(ext) -len(ext2)], OUTTAG, ext2, ext)

	scanner, file := utils.ReturnReader(filename, 0, false)
	defer file.Close()
	writer := utils.ReturnWriter(outfile, COMPRESSIONMODE, false)
	defer writer.Close()
	notCheckBarcode := true


	if barcodefilename != "" {
		notCheckBarcode = false
		refbarcodes = loadCellIDDict(barcodefilename)
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
			if refbarcodes[barcode] || notCheckBarcode {
				tocopy = true
				nbLine++
			}
		}

		if tocopy {
			if (isfour == 0) && (len(TAG) > 0) {
				buffer.WriteRune('@')
				buffer.WriteString(barcode)
				buffer.WriteString(TAG)
				buffer.WriteRune(':')
				buffer.WriteString(readHeader[1])
				buffer.WriteRune('\n')

			} else {
				buffer.WriteString(line)
				buffer.WriteRune('\n')
			}

			writer.Write(buffer.Bytes())
			buffer.Reset()
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

/*writeComplement ...*/
func writeComplement(filename string, complStrategy string) (nbLines int) {
	var buffer bytes.Buffer
	ext := path.Ext(filename)
	ext2 := path.Ext(filename[:len(filename) - len(ext)])

	outfile := fmt.Sprintf("%s.compl%s%s",
		filename[:len(filename) - len(ext) -len(ext2)], ext2, ext)

	scanner, file := utils.ReturnReader(filename, 0, false)
	defer file.Close()
	writer := utils.ReturnWriter(outfile, COMPRESSIONMODE, false)
	defer writer.Close()

	cycle := 0

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

				buffer.WriteRune('@')
				buffer.WriteString(split1)
				buffer.Write(returnComp(split2))
				buffer.WriteRune(':')
				buffer.WriteString(rest)
				buffer.WriteRune('\n')

			case "split_10_compl_first":
				split1 := barcode[:10]
				split2 := barcode[10:]

				buffer.WriteRune('@')
				buffer.Write(returnComp(split1))
				buffer.WriteString(split2)
				buffer.WriteRune(':')
				buffer.WriteString(rest)
				buffer.WriteRune('\n')

			default:
				panic(fmt.Sprintf("wrong strategy used!: %s\n", COMPLSTRATEGY))
			}

		default:
			buffer.WriteString(line)
			buffer.WriteRune('\n')

		}
		writer.Write(buffer.Bytes())
		buffer.Reset()

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
	var buffer bytes.Buffer

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
		buffer.WriteString(key)
		buffer.WriteString(SEP)
		buffer.WriteString(strconv.Itoa(value))
		buffer.WriteRune('\n')
		outfile.Write(buffer.Bytes())
		buffer.Reset()
	}

	if SORTLOGS {
		defer utils.SortLogfile(outfilename, SEP, "", IGNORESORTINGCATEGORY, IGNOREERROR)
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
	writer := utils.ReturnWriterForBzipfile("tmp.bz2")

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
		defer utils.SortLogfile(outfname, SEP, "", IGNORESORTINGCATEGORY, IGNOREERROR)
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

func loadCellIDDict(fname string) map[string]bool {
	f, err := os.Open(fname)
	check(err)
	defer f.Close()
	scanner := bufio.NewScanner(f)

	celliddict := make(map[string]bool)

	for scanner.Scan() {
		line := scanner.Text()

		celliddict[line] = true
	}
	return celliddict
}
