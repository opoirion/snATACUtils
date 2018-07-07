package ATACdemultiplexUtils

import "bufio"
import "io"
import "os"
import "log"
import "path"
import "github.com/dsnet/compress/bzip2"
import "os/exec"
import "strings"
import "compress/gzip"
import "fmt"

import (
	originalbzip2  "compress/bzip2"
	Cbzip2 "ATACdemultiplex/cbzip2"
)


const BUFFERSIZE = 1000000


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

/*ReturnWriterForBzipfile ... */
func ReturnWriterForBzipfile(fname string, compressionMode int) (io.WriteCloser) {
	outputFile, err := os.Create(fname)
	Check(err)
	bzipFile := Cbzip2.NewWriter(outputFile, compressionMode)
	Check(err)

	return bzipFile
}

/*ReturnWriter ... */
func ReturnWriter(fname string, compressionMode int, pureGo bool) (io.WriteCloser) {

	ext := path.Ext(fname)
	var bzipFile io.WriteCloser

	switch ext {
	case  ".bz2":
		switch pureGo {
		case true:
			bzipFile = ReturnWriterForBzipfilePureGo(fname)
		default:
			bzipFile = ReturnWriterForBzipfile(fname, compressionMode)
		}
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

/*ReturnWriterForBzipfilePureGo ... */
func ReturnWriterForBzipfilePureGo(fname string) (*bzip2.Writer) {
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

	readerBzip, fileOpen := returnCbzipReader(fname)

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
			readerBzip, fileOpen = returnCbzipReader(fname)
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


/*returnCbzipReader ... */
func returnCbzipReader(fname string) (io.ReadCloser, *os.File) {
	fileOpen, err := os.OpenFile(fname, 0, 0)

	if err != nil {
		log.Fatal(err)
	}

	readerOs := bufio.NewReader(fileOpen)
	readerBzip := Cbzip2.NewReader(readerOs)

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
