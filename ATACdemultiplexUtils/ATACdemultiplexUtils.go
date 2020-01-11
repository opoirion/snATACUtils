package atacdemultiplexutils

import (
	"bufio"
	"io"
	"os"
	"bytes"
	"path"
	"github.com/dsnet/compress/bzip2"
	"os/exec"
	"strings"
	gzip "github.com/klauspost/pgzip"
	"fmt"
	"sort"
	"strconv"
	"time"
	originalbzip2  "compress/bzip2"
)


/*Filename type used to check if files exists */
type Filename string

/*Set ... */
func (i *Filename) Set(filename string) error {
	if _, err := os.Stat(filename); err != nil {
		panic(fmt.Sprintf("!!!!Error %s with file: %s\n", err, filename))
	}

	*i =  Filename(filename)
	return nil
}


func (i *Filename) String() string {
	return string(*i)
}


/*ReturnReader Return reader for file */
func (i *Filename) ReturnReader(startingLine int) (*bufio.Scanner, *os.File) {
	return ReturnReader(string(*i), startingLine)
}

type closer interface {
	Close() error
}

/*BUFFERSIZE ... */
const BUFFERSIZE = 1000000

/*Pair ...*/
type Pair struct {
  Key string
  Value int
}

/*PairList ...*/
type PairList []Pair

// Len
func (p PairList) Len() int { return len(p) }
func (p PairList) Less(i, j int) bool { return p[i].Value < p[j].Value }
func (p PairList) Swap(i, j int){ p[i], p[j] = p[j], p[i] }

/*ArrayFlags ... */
type ArrayFlags []string

/*String ... */
func (i *ArrayFlags) String() string {
	var str string;
	for _, s := range (*i) {
		str = str + "\t"+  s
	}

	return str
}

/*Set ... */
func (i *ArrayFlags) Set(value string) error {
	*i = append(*i, value)
	return nil
}

/*RankByWordCountAndDeleteOldMap ...*/
func RankByWordCountAndDeleteOldMap(wordFrequencies * map[string]int) PairList{
	pl := make(PairList, len(*wordFrequencies))
	i := 0
	for k, v := range *wordFrequencies {
		pl[i] = Pair{k, v}
		i++
		delete(*wordFrequencies, k)
	}

	sort.Sort(sort.Reverse(pl))
	return pl
}


/*Check ... */
func Check(err error) {
	if err != nil {
		panic(err)
	}
}

/*CloseFile close file checking error */
func CloseFile(file closer) {
	err := file.Close()
	Check(err)
}

/*ExceCmd ... */
func ExceCmd(cmd string) {
	_, err := exec.Command("sh", "-c", cmd).Output()
	if err != nil {
		fmt.Printf("Error with cmd: %s\n", cmd)
	}

	Check(err)
}

/*ReturnWriter ... */
func ReturnWriter(fname string) (io.WriteCloser) {

	ext := path.Ext(fname)
	var bzipFile io.WriteCloser
	var err error

	switch ext {
	case  ".bz2":
		bzipFile = ReturnWriterForBzipfile(fname)

	case ".gz":
		bzipFile = ReturnWriterForGzipFile(fname)
	default:
		bzipFile, err = os.Create(fname)
		Check(err)
	}

	return bzipFile
}


/*ReturnWriterForGzipFile ... */
func ReturnWriterForGzipFile(fname string) (io.WriteCloser) {
	outputFile, err := os.Create(fname)
	Check(err)
	bzipFile := gzip.NewWriter(outputFile)
	Check(err)

	return bzipFile
}

/*ReturnWriterForBzipfile ... */
func ReturnWriterForBzipfile(fname string) (*bzip2.Writer) {
	outputFile, err := os.Create(fname)
	Check(err)
	bzipFile, err := bzip2.NewWriter(outputFile, new(bzip2.WriterConfig))
	Check(err)

	return bzipFile
}

/*ReturnReaderForBzipfileOld ... */
func ReturnReaderForBzipfileOld(fname string, seekPos int) (*bufio.Scanner, *os.File) {
	fileOpen, err := os.OpenFile(fname, 0, 0)
	Check(err)

	config := new(bzip2.ReaderConfig)

	readerOs := bufio.NewReader(fileOpen)
	readerBzip, err := bzip2.NewReader(readerOs, config)
	Check(err)
	readerOs2 := bufio.NewReader(readerBzip)
	bzipScanner := bufio.NewScanner(readerOs2)

	return bzipScanner, fileOpen
}

/*ReturnReader ... */
func ReturnReader(fname string, startingLine int) (*bufio.Scanner, *os.File) {
	ext := path.Ext(fname)
	var bzipScanner * bufio.Scanner
	var fileOpen * os.File
	var err error

	switch ext {
	case ".bz2":
		bzipScanner, fileOpen = ReturnReaderForBzipfile(fname, startingLine)

	case ".gz":
		bzipScanner, fileOpen = ReturnReaderForGzipfile(fname, startingLine)
	default:
		fileOpen, err = os.Open(fname)
		Check(err)
		bzipScanner = bufio.NewScanner(fileOpen)
	}

	return bzipScanner, fileOpen

}

/*ReturnReaderForGzipfile ... */
func ReturnReaderForGzipfile(fname string, startingLine int) (*bufio.Scanner, *os.File) {
	fileOpen, err := os.OpenFile(fname, 0, 0)
	Check(err)

	readerOs := bufio.NewReader(fileOpen)
	readerBzip, _ := gzip.NewReader(readerOs)
	bzipScanner := bufio.NewScanner(readerBzip)

	if startingLine > 0 {
		scanUntilStartingLine(bzipScanner, startingLine)
	}

	return bzipScanner, fileOpen
}

/*ReturnReaderForBzipfilePureGo ... */
func ReturnReaderForBzipfilePureGo(fname string, startingLine int) (*bufio.Scanner, *os.File) {
	var bzipScanner * bufio.Scanner
	buffer := make([]byte, BUFFERSIZE)

	fileOpen, err := os.OpenFile(fname, 0, 0)
	Check(err)

	readerOs := bufio.NewReader(fileOpen)
	readerBzip := originalbzip2.NewReader(readerOs)

	if startingLine == 0 {
		bzipScanner = bufio.NewScanner(readerBzip)
		return bzipScanner, fileOpen
	}

	_, err = readerBzip.Read(buffer)
	Check(err)

	nbLines := strings.Count(string(buffer), "\n")
	currentLine := nbLines

loop:
	for {
		switch {
		case nbLines > startingLine:
			fileOpen, _ := os.OpenFile(fname, 0, 0)
			readerOs := bufio.NewReader(fileOpen)
			readerBzip := originalbzip2.NewReader(readerOs)
			bzipScanner = bufio.NewScanner(readerBzip)
			scanUntilStartingLine(bzipScanner, startingLine)
			break loop

		default:
			bzipScanner = bufio.NewScanner(readerBzip)
			scanUntilStartingLine(bzipScanner, startingLine - currentLine)
			break loop
		}
	}
	return bzipScanner, fileOpen
}


/*ReturnReaderForBzipfile ... */
func ReturnReaderForBzipfile(fname string, startingLine int) (*bufio.Scanner, *os.File) {
	var bzipScanner * bufio.Scanner
	var err error

	buffer := make([]byte, BUFFERSIZE)

	readerBzip, fileOpen := returnBzipReader(fname)

	if startingLine == 0 {
		bzipScanner = bufio.NewScanner(readerBzip)
		return bzipScanner, fileOpen
	}

	_, err = readerBzip.Read(buffer)
	Check(err)

	nbLines := strings.Count(string(buffer), "\n")
	currentLine := nbLines

loop:
	for {
		switch {
		case nbLines > startingLine:
			readerBzip, fileOpen = returnBzipReader(fname)
			bzipScanner = bufio.NewScanner(readerBzip)
			scanUntilStartingLine(bzipScanner, startingLine)
			break loop

		default:
			bzipScanner = bufio.NewScanner(readerBzip)
			scanUntilStartingLine(bzipScanner, startingLine - currentLine)
			break loop
		}
	}
	return bzipScanner, fileOpen
}

/*scanUntilStartingLine ... */
func scanUntilStartingLine(scanner * bufio.Scanner, nbLine int) {
	var ok bool
	for i := 0;i < nbLine; i++ {
		ok = scanner.Scan()

		if !ok {
			break
		}
	}

}


/*returnBzipReader ... */
func returnBzipReader(fname string) (io.Reader, *os.File) {
	fileOpen, err := os.OpenFile(fname, 0, 0)
	Check(err)

	readerOs := bufio.NewReader(fileOpen)
	readerBzip := originalbzip2.NewReader(readerOs)

	return readerBzip, fileOpen
}


/*SortLogfile sort a file according to key value
input:
    filename string,
    separator string,
    outfname string,
    ignoreSortingCategory bool,
    ignoreError bool
*/
func SortLogfile(filename Filename, separator string, outfname string,
	ignoreSortingCategory bool, ignoreError bool)  {
	scanner, file := filename.ReturnReader(0)
	defer CloseFile(file)
	var buffer bytes.Buffer
	var split []string
	var valueField, key string
	var value int

	ext := path.Ext(filename.String())

	if outfname == "" {
		outfname = fmt.Sprintf("%s_sorted%s", strings.TrimSuffix(filename.String(), ext), ext)
	}

	outfile, err := os.Create(outfname)
	Check(err)
	defer CloseFile(outfile)
	defer os.Rename(outfname, filename.String())

	Check(err)
	pl := PairList{}

	lineNb := -1
	buff := 0

	for scanner.Scan() {
		line := scanner.Text()
		lineNb++

		if len(line) == 0 || line[0] == '#' || line[0] == '\n'  {
			if ignoreSortingCategory{
				continue
			}
			buffer.WriteString(line)
			buffer.WriteRune('\n')

			outfile.Write(buffer.Bytes())
			buffer.Reset()

			if len(pl) == 0  {
				continue
			}
			sort.Slice(pl, func(i, j int) bool {
				return pl[i].Value > pl[j].Value
			})

			for _, el := range pl {
				buffer.WriteString(el.Key)
				buffer.WriteString(separator)
				buffer.WriteString(strconv.Itoa(el.Value))
				buffer.WriteRune('\n')

				outfile.Write(buffer.Bytes())
				buffer.Reset()

				buff++
				if buff > 100000{
					outfile.Sync()
					buff = 0
				}
			}

			pl = PairList{}
			continue
		}

		split = strings.Split(line, separator)
		valueField = split[len(split)-1]

		if len(split) > 2 {
			key = strings.Join(split[:len(split)-1], separator)
		} else {
			key = split[0]
		}
		value, err = strconv.Atoi(valueField)

		if err != nil {
			fmt.Printf("value field %s from: %s at line nb %d not conform!\n",
				valueField, line, lineNb)

			if !ignoreError {
				Check(err)
			}
		}


		pl = append(pl, Pair{Key:key, Value:value})
	}

	if len(pl) == 0 {
		return
	}

	sort.Slice(pl, func(i, j int) bool {
		return pl[i].Value > pl[j].Value
	})

	for _, el := range pl {
		outfile.WriteString(fmt.Sprintf("%s%s%d\n", el.Key, separator, el.Value))
		buff++
		if buff > 100000{
			err = outfile.Sync()
			Check(err)
			buff = 0
		}
	}
}

/*LoadCellIDDict create cell ID bool dict */
func LoadCellIDDict(fname string) map[string]bool {
	f, err := os.Open(fname)
	Check(err)
	defer CloseFile(f)
	scanner := bufio.NewScanner(f)
	var cellID string

	celliddict := make(map[string]bool)

	for scanner.Scan() {
		cellID = scanner.Text()
		cellID = strings.ReplaceAll(cellID, " ", "\t")
		cellID = strings.Split(cellID, "\t")[0]

		celliddict[cellID] = true
	}
	return celliddict
}


/*LoadCellDictsToIndex create cell ID index dict */
func LoadCellDictsToIndex(fname Filename) map[string]int {

	scanner, f := fname.ReturnReader(0)
	defer CloseFile(f)

	var cellID string
	var index int

	celliddict := make(map[string]int)

	for scanner.Scan() {
		cellID = scanner.Text()
		cellID = strings.ReplaceAll(cellID, " ", "\t")
		cellID = strings.Split(cellID, "\t")[0]

		celliddict[cellID] = index
		index++
	}

	return celliddict
}


/*LoadCellDictsFromBedFileToIndex create cell ID <-> index dict using a single cell bed file storing unique CELL ID
Ex:
chr1    15555    15522    CELLID1
chr1    25555    25522    CELLID1
chr1    35555    35522    CELLID1
chr1    15555    15522    CELLID2
chr1    25555    25522    CELLID2

will give map(CELLID1:0, CELLID2:1)
  */
func LoadCellDictsFromBedFileToIndex(fname Filename) map[string]int {

	scanner, f := fname.ReturnReader(0)
	defer CloseFile(f)

	var cellID, line string
	var index int
	var isInside bool

	celliddict := make(map[string]int)

	tStart := time.Now()

	for scanner.Scan() {
		line = scanner.Text()
		cellID = strings.Split(line, "\t")[3]

		if _, isInside = celliddict[cellID];!isInside {
			celliddict[cellID] = index
			index++

		}
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Loading Cell Dicts From Bed File To Index done in time: %f s \n", tDiff.Seconds())

	return celliddict
}


/*CountNbLines count nb lines in a file*/
func CountNbLines(filename string) int {
	reader, file := ReturnReader(filename, 0)
	defer CloseFile(file)

	nbLines := 0

	tStart := time.Now()

	for reader.Scan() {
		nbLines++
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Count nb lines done in time: %f s \n", tDiff.Seconds())

	return nbLines
}
