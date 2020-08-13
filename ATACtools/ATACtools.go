/* Suite of functions dedicated to pre/post process files related to snATAC pipeline */

package main

import(
	"fmt"
	"flag"
	"bufio"
	"os"
	"log"
	"time"
	"strings"
	"path"
	"sync"
	"strconv"
	"bytes"
	// "github.com/dsnet/compress/bzip2"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
)

/*FILENAME ...*/
var FILENAME utils.Filename

/*MAX ...*/
var MAX int

/*PRINTLASTLINE ...*/
var PRINTLASTLINE bool

/*PRINTLASTLINES ...*/
var PRINTLASTLINES int

/*GOTOLINE ...*/
var GOTOLINE int

/*PATTERN ...*/
var PATTERN string

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

/*SCAN ...*/
var SCAN bool

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

/*CLEAN clean files with unwanted lines*/
var CLEAN bool

/*COUNT count patterns in file*/
var COUNT bool

/*CLEANPATTERN pattern used to clean files with unwanted lines*/
var CLEANPATTERN string

/*MAXSCANTOKENSIZE int*/
var MAXSCANTOKENSIZE int

/*THREADNB int*/
var THREADNB int

func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `
#################### Suite of functions dedicated to pre/post process files related to snATAC pipeline ########################

-bed_to_cicero: format a bed to cicero input (ex: chr1\t1215523\t1216200\tcellID -> chr1_1215523_121620,tcellID,1)
USAGE: ATACtools -bed_to_cicero -filename <bedfile> (-filenames <bedfile2> -filenames <bedfile3> ...  -ignoreerror)

-create_ref_fastq: Create a ref FASTQ file using a reference barcode list
USAGE: ATACtools -create_ref_bed -filename <fname> (-ref_barcode_list <fname> -tag <string> -output <string>)

-merge: Merge input log files together
USAGE: ATACtools -merge -filenames <fname1> -filenames <fname2> -filenames <fname3> ...  (-sortfile -delimiter "<string>" -ignoreerror -ignore_sorting_category)

-sortfile: Sort key -> values (i.e.: <key><SEP><value>) file
USAGE: ATACtools -sortfile -filename <fname> (-delimiter <string> -ignoreerror -ignore_sorting_category)

-write_compl: Write the barcode complement of a fastq files
USAGE: ATACtools -write_compl <fastq_file> (-compl_strategy <"split_10_compl_second"/"split_10_compl_first"> -tag <string>)

-scan: Scan a file and determine the number of line
USAGE ATACtools -scan -filename <string> (-printlastline -printlastlines <int> -pattern <string> -gotoline <int>)

-create_barcode_dict: Create a barcode key / value count file
USAGE: ATACtools -create_barcode_dict -filename <fname> (-sortfile -delimiter <string>)

-clean: clean file from unwanted lines
USAGE: ATACtools -clean -filename <fname> -output filename -clean_pattern "\n"

-count: Count string in file (useful for large file)
USAGE: ATACtools -count -filename (-pattern <string> -threads <int>)
`)
		 flag.PrintDefaults()
	}


	flag.Var(&FILENAME, "filename", "name of the file(s) multiple files should be into \"")
	flag.Var(&FILENAMES, "filenames", "name of the files to use")
	flag.IntVar(&MAX, "max_nb_lines", 0, "max number of lines")
	flag.BoolVar(&CICEROPROCESSING, "bed_to_cicero", false,
		`format a bed to cicero input (ex: chr1\t1215523\t1216200\tcellID -> chr1_1215523_121620,tcellID,1)`)
	flag.BoolVar(&CREATEREFFASTQ, "create_ref_fastq", false, `create a ref FASTQ file using a reference barcode list`)
	flag.BoolVar(&CREATEREFBEDFILE, "create_ref_bed", false, `create a ref bed file using a reference barcode list`)
	flag.StringVar(&REFBARCODELIST, "ref_barcode_list", "", "file containing the reference barcodes (one per line)")
	flag.BoolVar(&MERGE, "merge", false, `merge input log files together`)
	flag.IntVar(&GOTOLINE, "gotoline", 0, "go to line")
	flag.BoolVar(&SORTLOGS, "sortfile", false,
		`sort files (<key><SEP><value> file`)
	flag.BoolVar(&PRINTLASTLINE, "printlastline", false, "print last line")
	flag.IntVar(&PRINTLASTLINES, "printlastlines", 0, "print last n lines")
	flag.BoolVar(&IGNOREERROR, "ignoreerror", false, "ignore error and continue")
	flag.BoolVar(&IGNORESORTINGCATEGORY, "ignore_sorting_category", false, "ignore file cateogry (identified by #) when sorting")
	flag.StringVar(&PATTERN, "pattern", "", "search specific string in file")
	flag.StringVar(&OUTFILE, "output", "", "file name of the output")
	flag.StringVar(&OUTTAG, "output_tag", "reference", "particule to annotate the output file name")
	flag.BoolVar(&WRITECOMPL, "write_compl", false, `write the barcode complement of a fastq files`)
	flag.BoolVar(&SCAN, "scan", false, `scan a file and determine the number of line`)
	flag.BoolVar(&COUNT, "count", false, `Count patterns in file`)
	flag.StringVar(&SEP, "delimiter", "\t", "delimiter used to split and sort the log file (default \t)")
	flag.StringVar(&TAG, "tag", "", "tag used when creating a reference fastq file to tag all the reads (default \"\")")
	flag.StringVar(&COMPLSTRATEGY, "compl_strategy", "split_10_compl_second", `Strategy to use when writing the complement of a fastq file (default split_10_compl_second: split after 10 bases and complementary only second)`)
	flag.BoolVar(&CREATEBARCODEDICT, "create_barcode_dict", false, `create a barcode key / value count file`)
	flag.BoolVar(&CLEAN, "clean", false, `clean files with unwanted lines`)
	flag.StringVar(&CLEANPATTERN, "clean_pattern", "\n", "pattern used to clean files with unwanted lines")
	flag.IntVar(&MAXSCANTOKENSIZE, "max_scan_size", 0, "MaxScanTokenSize variable that defines the length of a line that a buffer can read (set if differnt than 0)")
	flag.IntVar(&THREADNB, "threads", 8, "threads concurrency for specific usage")
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
			go extractFASTQreadsPerBarcodes(utils.Filename(filename), REFBARCODELIST, &waiting)
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
	case SCAN:
		nbLines = countLine(FILENAME)
	case CLEAN:
		cleanFile(FILENAME)
	case COUNT:
		countInFile(FILENAME)
	default:
		log.Fatal("Error at least one processing option (scan/bed_to_cicero/create_ref_fastq/create_ref_bed/merge/create_barcode_list) should be used!")
	}

	tDiff := time.Since(tStart)

	if nbLines > 0 || SCAN {
		fmt.Printf("number of lines: %d\n", nbLines)
	}

	fmt.Printf("time: %f s\n", tDiff.Seconds())

}

func countInFile(filename utils.Filename) {
	_, file := filename.ReturnReader(0)
	defer utils.CloseFile(file)

	buffersize := 10000000
	buffers := make([][]byte, THREADNB)
	availableBuffers := make(chan int, THREADNB)
	var reader utils.Reader

	countMap := make([]int, THREADNB)
	count := 0

	if PATTERN == "" {
		PATTERN = "\n"
	}

	pattern := []byte(PATTERN)
	waiting := sync.WaitGroup{}

	var err error
	var index, nbbytes int

	for i := 0; i < THREADNB;i++ {
		availableBuffers <- i
		buffers[i] = make([]byte, buffersize)
	}

	reader = utils.ReturnFileReader(filename.String())
	var remakeBuffer int

	for {
		index = <- availableBuffers
		nbbytes, err = reader.Read(buffers[index])

		if err != nil {
			break
		}

		if nbbytes < buffersize {
			buffers[index] = buffers[index][:nbbytes]
			remakeBuffer = buffersize
		} else {
			remakeBuffer = 0
		}

		waiting.Add(1)

		go countAndUpdate(
			pattern,
			&buffers[index],
			index,
			&countMap,
			remakeBuffer,
			&waiting,
			availableBuffers)
	}

	waiting.Wait()

	for _, c := range countMap {
		count += c
	}
	fmt.Printf("Scanning file: %s done\n (status: %s) \n",
		filename.String(), err )
	fmt.Printf("Count of pattern (%q): %d \n", PATTERN, count )
}

func countAndUpdate(pattern []byte,
	buffers * []byte,
	index int,
	countMap * []int,
	remakeBuffer int,
	waiting * sync.WaitGroup,
	availableBuffers chan int) {
	defer waiting.Done()
	countLocal := bytes.Count((*buffers), pattern)

	(*countMap)[index] += countLocal

	if remakeBuffer > 0 {
		(*buffers) = make([]byte, remakeBuffer)
	}

	availableBuffers <- index
}

func cleanFile(filename utils.Filename) {
	var line string
	var err error
	scanner, file := filename.ReturnReader(0)
	defer utils.CloseFile(file)

	if OUTFILE == "" {
		fsplit := strings.SplitN(filename.String(), ".", 2)
		OUTFILE = fmt.Sprintf("%s.clean.%s", fsplit[0], fsplit[1])
	}

	writer := utils.ReturnWriter(OUTFILE)
	defer utils.CloseFile(writer)
	var buffer bytes.Buffer
	count := 0
	countrm := 0

	for scanner.Scan() {
		line = scanner.Text()

		if line == CLEANPATTERN {
			countrm++
			continue
		}

		buffer.WriteString(line)
		buffer.WriteRune('\n')
		count++

		if count >= 50000 {
			_, err = writer.Write(buffer.Bytes())
			utils.Check(err)
			count = 0
			buffer.Reset()
		}
	}

	_, err = writer.Write(buffer.Bytes())
	utils.Check(err)

	fmt.Printf("%d lines removed. output file: %s written!\n", countrm, OUTFILE)
}


/*bedFilestoCiceroInput ...*/
func bedFilestoCiceroInput(filenames []string, outfile string) {
	var ext string
	var finalOutfile string

	if FILENAME != "" {
		filenames = append(filenames, FILENAME.String())
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

	tDiff := time.Since(tStart)
	fmt.Printf("file:%s created in %f s\n", finalOutfile, tDiff.Seconds())
}

/*bedFiletoCiceroInputOnethread ...*/
func bedFiletoCiceroInputOnethread(filename string, outfile string, barcodefilename string,
	waiting * sync.WaitGroup, writeHeader bool) {
	defer waiting.Done()
	var buffer bytes.Buffer
	var split = make([]string, 4)
	var nbLine int
	var barcodeIndex map[string]bool
	var checkBarcode bool

	if barcodefilename != "" {
		checkBarcode = true
		barcodeIndex = utils.LoadCellIDDict(barcodefilename)
	}

	scanner, file := utils.ReturnReader(filename, 0)
	defer file.Close()
	writer := utils.ReturnWriter(outfile)
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
func extractBEDreadsPerBarcodes(filename utils.Filename, barcodefilename string) {

	var barcode string
	var refbarcodes map[string]bool
	var split = make([]string, 4)
	var notCheckBarcode = true
	var buffer bytes.Buffer
	var outfile string

	ext := path.Ext(filename.String())
	ext2 := path.Ext(filename.String()[:len(filename) - len(ext)])

	if OUTFILE != "" {
		outfile = OUTFILE
	} else {
		outfile = fmt.Sprintf("%s.%s%s%s",
			filename[:len(filename) - len(ext) -len(ext2)], OUTTAG, ext2, ext)
	}

	scanner, file := filename.ReturnReader(0)
	defer file.Close()
	writer := utils.ReturnWriter(outfile)
	defer writer.Close()

	if barcodefilename != "" {
		notCheckBarcode = false
		refbarcodes = utils.LoadCellIDDict(barcodefilename)
	}

	nbLine := 0
	buffCount := 0

	for scanner.Scan() {
		line := scanner.Text()
		nbLine++

		split = strings.Split(line, "\t")
		barcode = split[3]

		if notCheckBarcode || refbarcodes[barcode] {

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

			buffCount++

			if buffCount > 100000 {
				writer.Write(buffer.Bytes())
				buffer.Reset()
				buffCount = 0
			}
		}
	}

	writer.Write(buffer.Bytes())
	buffer.Reset()

	fmt.Printf("nb barcodes extracted: %d\n", nbLine)
	fmt.Printf("file created: %s\n", outfile)
}

/*extractFASTQreadsPerBarcodes create a new FASTQ file using a barcode name */
func extractFASTQreadsPerBarcodes(filename utils.Filename, barcodefilename string, waiting *sync.WaitGroup) {

	var barcode string
	var refbarcodes = make(map[string]bool)
	var readHeader = make([]string, 2)
	defer waiting.Done()
	var buffer bytes.Buffer
	var outfile string

	ext := path.Ext(filename.String())
	ext2 := path.Ext(filename.String()[:len(filename) - len(ext)])

	if OUTFILE != "" {
		outfile = OUTFILE
	} else {

		outfile = fmt.Sprintf("%s.%s%s%s",
		filename[:len(filename) - len(ext) -len(ext2)], OUTTAG, ext2, ext)
	}

	scanner, file := filename.ReturnReader(0)
	defer file.Close()
	writer := utils.ReturnWriter(outfile)
	defer writer.Close()
	notCheckBarcode := true


	if barcodefilename != "" {
		notCheckBarcode = false
		refbarcodes = utils.LoadCellIDDict(barcodefilename)
	}

	isfour := 0
	nbLine := 0
	tocopy := false
	buffCount := 0

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

			buffCount++

			if buffCount > 100000 {
				writer.Write(buffer.Bytes())
				buffer.Reset()
				buffCount = 0
			}
		}

		isfour++

		if isfour == 4 {
			tocopy = false
			isfour = 0
		}
	}

	writer.Write(buffer.Bytes())
	buffer.Reset()

	fmt.Printf("nb barcodes extracted: %d\n", nbLine)
	fmt.Printf("file created: %s\n", outfile)
}

/*writeComplement ...*/
func writeComplement(filename utils.Filename, complStrategy string) (nbLines int) {
	var buffer bytes.Buffer
	ext := path.Ext(filename.String())
	ext2 := path.Ext(filename.String()[:len(filename) - len(ext)])
	var outfile string

	if OUTFILE != "" {
		outfile = OUTFILE
	} else {
		outfile = fmt.Sprintf("%s.compl%s%s",
			filename[:len(filename) - len(ext) -len(ext2)], ext2, ext)
	}

	scanner, file := filename.ReturnReader(0)
	defer file.Close()
	writer := utils.ReturnWriter(outfile)
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
		defer createIndexCountFile(utils.Filename(outfile))
	}

	return nbLines
}

func createIndexCountFile(filename utils.Filename) (nbLines int) {

	ext := path.Ext(filename.String())
	ext2 := path.Ext(filename.String()[:len(filename) - len(ext)])
	var outfilename string

	if OUTFILE != "" {
		outfilename = OUTFILE
	} else {
		outfilename = fmt.Sprintf("%s.barcodeCounts",
			filename[:len(filename) - len(ext) -len(ext2)])
	}

	scanner, file := filename.ReturnReader(0)
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
		defer utils.SortLogfile(utils.Filename(outfilename), SEP, "", IGNORESORTINGCATEGORY, IGNOREERROR)
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

func countLine(filename utils.Filename) int {
	scanner, file := filename.ReturnReader(0)
	defer file.Close()
	return processScanner(scanner)
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
		utils.Check(err)

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
		defer utils.SortLogfile(
			utils.Filename(outfname), SEP, "", IGNORESORTINGCATEGORY, IGNOREERROR)
	}
}

func processScanner(scanner * bufio.Scanner) (nbLines int) {
	nbLines = 0


	if MAXSCANTOKENSIZE > 0 {
		scanner.Buffer([]byte{}, MAXSCANTOKENSIZE * MAXSCANTOKENSIZE)
	}


	for scanner.Scan() {
		nbLines++

		if (GOTOLINE > 0) && (PRINTLASTLINES > 0) && (GOTOLINE - nbLines < PRINTLASTLINES)  {
			line := scanner.Text()
			fmt.Printf("last lines: %s\n", line)
		}

		if GOTOLINE > 0 && nbLines > GOTOLINE {
			break
		}

		if PATTERN != "" {
			line := scanner.Text()
			if strings.Contains(line, PATTERN) {
				fmt.Printf("line nb: %d\t %s found in %s\n", nbLines, PATTERN, line)
			}
		}
	}

	if PRINTLASTLINE {
		line := scanner.Text()
		fmt.Printf("last line: %s\n", line)
	}

	utils.Check(scanner.Err())

	return nbLines

}


func check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
