package ATACdemultiplexUtils

import (
	"bufio"
	"io"
	"os"
	"log"
	"bytes"
	"path"
	"github.com/dsnet/compress/bzip2"
	"os/exec"
	"strings"
	gzip "github.com/klauspost/pgzip"
	"fmt"
	"sort"
	"strconv"
	originalbzip2  "compress/bzip2"
)


const BUFFERSIZE = 1000000

/*Pair ...*/
type Pair struct {
  Key string
  Value int
}

/*PairList ...*/
type PairList []Pair

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
		log.Fatal(err)
	}
}


/*ExceCmd ... */
func ExceCmd(cmd string) {
	_, err := exec.Command("sh", "-c", cmd).Output()
	Check(err)
}

/*ReturnWriter ... */
func ReturnWriter(fname string, compressionMode int, pureGo bool) (io.WriteCloser) {

	ext := path.Ext(fname)
	var bzipFile io.WriteCloser

	switch ext {
	case  ".bz2":
		bzipFile = ReturnWriterForBzipfile(fname)

	case ".gz":
		bzipFile = ReturnWriterForGzipFile(fname)
	default:
		panic(fmt.Sprintf("%s does not have either bzip2 (bz) or gzip (gz) extension!",
			fname))
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

	if err != nil {
		log.Fatal(err)
	}
	config := new(bzip2.ReaderConfig)

	readerOs := bufio.NewReader(fileOpen)
	readerBzip, err := bzip2.NewReader(readerOs, config)
	Check(err)
	readerOs2 := bufio.NewReader(readerBzip)
	bzipScanner := bufio.NewScanner(readerOs2)

	return bzipScanner, fileOpen
}

/*ReturnReader ... */
func ReturnReader(fname string, startingLine int, pureGo bool) (*bufio.Scanner, *os.File) {
	ext := path.Ext(fname)
	var bzipScanner * bufio.Scanner
	var fileOpen * os.File

	switch ext {
	case ".bz2":
		switch pureGo {
		case true:
			bzipScanner, fileOpen = ReturnReaderForBzipfilePureGo(fname, startingLine)
		default:
			bzipScanner, fileOpen = ReturnReaderForBzipfile(fname, startingLine)
		}
	case ".gz":
		bzipScanner, fileOpen = ReturnReaderForGzipfile(fname, startingLine)
	default:
		panic(fmt.Sprintf("%s does not have either bzip2 (bz) or gzip (gz) extension!",
			fname))
	}

	return bzipScanner, fileOpen

}

/*ReturnReaderForGzipfile ... */
func ReturnReaderForGzipfile(fname string, startingLine int) (*bufio.Scanner, *os.File) {
	fileOpen, err := os.OpenFile(fname, 0, 0)

	if err != nil {
		log.Fatal(err)
	}

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
	buffer := make([]byte, BUFFERSIZE, BUFFERSIZE)

	fileOpen, err := os.OpenFile(fname, 0, 0)

	if err != nil {
		log.Fatal(err)
	}

	readerOs := bufio.NewReader(fileOpen)
	readerBzip := originalbzip2.NewReader(readerOs)

	if startingLine == 0 {
		bzipScanner = bufio.NewScanner(readerBzip)
		return bzipScanner, fileOpen
	}

	readerBzip.Read(buffer)

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

	buffer := make([]byte, BUFFERSIZE, BUFFERSIZE)

	readerBzip, fileOpen := returnBzipReader(fname)

	if startingLine == 0 {
		bzipScanner = bufio.NewScanner(readerBzip)
		return bzipScanner, fileOpen
	}

	readerBzip.Read(buffer)

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

	if err != nil {
		log.Fatal(err)
	}

	readerOs := bufio.NewReader(fileOpen)
	readerBzip := originalbzip2.NewReader(readerOs)

	return readerBzip, fileOpen
}


/*seekFile ... */
func seekFile(reader * io.ReadCloser, pos int) {
	currentPos := 0
	buff := make([]byte, BUFFERSIZE)

loop:
	for {
		if pos - currentPos < BUFFERSIZE {
			buff := make([]byte, pos - currentPos)
			(*reader).Read(buff)
			break loop
		}

		_, err := (*reader).Read(buff)
		Check(err)
		currentPos += BUFFERSIZE
	}
}

/*SortLogfile sort a file according to key value
input:
    filename string,
    separator string,
    outfname string,
    ignoreSortingCategory bool,
    ignoreError bool
*/
func SortLogfile(filename string, separator string, outfname string,
	ignoreSortingCategory bool, ignoreError bool)  {
	file, err := os.Open(filename)
	defer file.Close()
	var buffer bytes.Buffer
	var split []string
	var valueField, key string
	var value int
	check(err)
	scanner := bufio.NewScanner(file)

	ext := path.Ext(filename)

	if outfname == "" {
		outfname = fmt.Sprintf("%s_sorted%s", strings.TrimSuffix(filename, ext), ext)
	}

	outfile, err := os.Create(outfname)
	defer outfile.Close()

	defer os.Rename(outfname, filename)

	check(err)
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
				check(err)
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
			outfile.Sync()
			buff = 0
		}
	}
}

func check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}

func LoadCellIDDict(fname string) map[string]bool {
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
