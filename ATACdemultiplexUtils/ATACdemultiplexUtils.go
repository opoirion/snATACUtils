package ATACdemultiplexUtils

import "bufio"
import "io"
import "os"
import "log"
import "github.com/dsnet/compress/bzip2"
import "os/exec"
import "strings"

import (
	originalbzip2  "compress/bzip2"
	Cbzip2 "ATACdemultiplex/cbzip2"
)


const BUFFER_SIZE = 1000000


func Check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}


func ExceCmd(cmd string) {
	_, err := exec.Command("sh", "-c", cmd).Output()
	Check(err)
}

func ReturnWriterForBzipfile(fname string, compressionMode int) (io.WriteCloser) {
	output_file, err := os.Create(fname)
	Check(err)
	bzip_file := Cbzip2.NewWriter(output_file, compressionMode)
	Check(err)

	return bzip_file
}

func ReturnWriterForBzipfilePureGo(fname string) (*bzip2.Writer) {
	output_file, err := os.Create(fname)
	Check(err)
	bzip_file, err := bzip2.NewWriter(output_file, new(bzip2.WriterConfig))
	Check(err)

	return bzip_file
}

func ReturnReaderForBzipfileOld(fname string, seekPos int) (*bufio.Scanner, *os.File) {
	file_open, err := os.OpenFile(fname, 0, 0)

	if err != nil {
		log.Fatal(err)
	}
	config := new(bzip2.ReaderConfig)

	reader_os := bufio.NewReader(file_open)
	reader_bzip, err := bzip2.NewReader(reader_os, config)
	Check(err)
	reader_os2 := bufio.NewReader(reader_bzip)
	bzip_scanner := bufio.NewScanner(reader_os2)

	return bzip_scanner, file_open
}


func ReturnReaderForBzipfilePureGoOld(fname string, seekPos int) (*bufio.Scanner, *os.File) {
	file_open, err := os.OpenFile(fname, 0, 0)

	if err != nil {
		log.Fatal(err)
	}

	reader_os := bufio.NewReader(file_open)
	reader_bzip := originalbzip2.NewReader(reader_os)

	// if seekPos > 0 {
	// 	seekFile(reader_bzip, seekPos)
	// }
	// reader_os2 := bufio.NewReader(reader_bzip)
	bzip_scanner := bufio.NewScanner(reader_bzip)

	return bzip_scanner, file_open
}

func ReturnReaderForBzipfilePureGo(fname string, startingLine int) (*bufio.Scanner, *os.File) {
	var bzip_scanner * bufio.Scanner
	buffer := make([]byte, BUFFER_SIZE, BUFFER_SIZE)

	file_open, err := os.OpenFile(fname, 0, 0)

	if err != nil {
		log.Fatal(err)
	}

	reader_os := bufio.NewReader(file_open)
	reader_bzip := originalbzip2.NewReader(reader_os)

	if startingLine == 0 {
		bzip_scanner = bufio.NewScanner(reader_bzip)
		return bzip_scanner, file_open
	}

	reader_bzip.Read(buffer)

	nbLines := strings.Count(string(buffer), "\n")
	currentLine := nbLines


	bzip_scanner = bufio.NewScanner(reader_bzip)
	scanUntilStartingLine(bzip_scanner, startingLine - currentLine)

	return bzip_scanner, file_open
}


func ReturnReaderForBzipfile(fname string, startingLine int) (*bufio.Scanner, *os.File) {
	var bzip_scanner * bufio.Scanner

	buffer := make([]byte, BUFFER_SIZE, BUFFER_SIZE)

	reader_bzip, file_open := returnCbzipReader(fname)

	if startingLine == 0 {
		bzip_scanner = bufio.NewScanner(reader_bzip)
		return bzip_scanner, file_open
	}

	reader_bzip.Read(buffer)

	nbLines := strings.Count(string(buffer), "\n")
	currentLine := nbLines

	bzip_scanner = bufio.NewScanner(reader_bzip)
	scanUntilStartingLine(bzip_scanner, startingLine - currentLine)

// loop:
// 	for {
// 		switch {
// 		case nbLines > startingLine:
// 			reader_bzip, file_open = returnCbzipReader(fname)
// 			bzip_scanner = bufio.NewScanner(reader_bzip)
// 			scanUntilStartingLine(bzip_scanner, startingLine)
// 			break loop

// 		case currentLine < startingLine - nbLines:
// 			reader_bzip.Read(buffer)
// 			nbLineIt := strings.Count(string(buffer), "\n")
// 			currentLine += nbLineIt

// 		default:
// 			bzip_scanner = bufio.NewScanner(reader_bzip)
// 			scanUntilStartingLine(bzip_scanner, startingLine - currentLine)
// 			break loop
// 		}
// 	}
	return bzip_scanner, file_open
}

func scanUntilStartingLine(scanner * bufio.Scanner, nbLine int) {
	for i := 0;i < nbLine; i++ {
		scanner.Scan()
	}

}


func returnCbzipReader(fname string) (io.ReadCloser, *os.File) {
	file_open, err := os.OpenFile(fname, 0, 0)

	if err != nil {
		log.Fatal(err)
	}

	reader_os := bufio.NewReader(file_open)
	reader_bzip := Cbzip2.NewReader(reader_os)

	return reader_bzip, file_open
}

func seekLine(fname, startingLine int)  {


}


func seekFile(reader * io.ReadCloser, pos int) {
	currentPos := 0
	buff := make([]byte, BUFFER_SIZE)

loop:
	for {
		if pos - currentPos < BUFFER_SIZE {
			buff := make([]byte, pos - currentPos)
			(*reader).Read(buff)
			break loop
		}

		_, err := (*reader).Read(buff)
		Check(err)
		currentPos += BUFFER_SIZE
	}
}
