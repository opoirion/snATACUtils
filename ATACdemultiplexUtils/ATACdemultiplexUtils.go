package ATACdemultiplexUtils

import "bufio"
import "io"
import "os"
import "log"
import "github.com/dsnet/compress/bzip2"
import "os/exec"

import (
	originalbzip2  "compress/bzip2"
	Cbzip2 "ATACdemultiplex/cbzip2"
)


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

func ReturnReaderForBzipfileOld(fname string, seekPos int64) (*bufio.Scanner, *os.File) {
	file_open, err := os.OpenFile(fname, 0, 0)
	file_open.Seek(seekPos, 0)

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


func ReturnReaderForBzipfilePureGo(fname string, seekPos int64) (*bufio.Scanner, *os.File) {
	file_open, err := os.OpenFile(fname, 0, 0)
	file_open.Seek(seekPos, 0)

	if err != nil {
		log.Fatal(err)
	}

	reader_os := bufio.NewReader(file_open)
	reader_bzip := originalbzip2.NewReader(reader_os)
	// reader_os2 := bufio.NewReader(reader_bzip)
	bzip_scanner := bufio.NewScanner(reader_bzip)

	return bzip_scanner, file_open
}

func ReturnReaderForBzipfile(fname string, seekPos int64) (*bufio.Scanner, *os.File) {
	file_open, err := os.OpenFile(fname, 0, 0)
	file_open.Seek(seekPos, 0)

	if err != nil {
		log.Fatal(err)
	}

	reader_os := bufio.NewReader(file_open)
	reader_bzip := Cbzip2.NewReader(reader_os)
	bzip_scanner := bufio.NewScanner(reader_bzip)


	return bzip_scanner, file_open
}
