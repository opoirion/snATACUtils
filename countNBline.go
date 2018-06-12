package main

import(
	"os";
	"bytes";
	"io"
)




var FILENAME string
var MAX int


func countLine(fname string) int {
	reader, err := os.Open(fname)
	check(err)

	buf := make([]byte, 32*1024)
	count := 0
	lineSep := []byte{'\n'}

	for {
		c, err := reader.Read(buf)
		count += bytes.Count(buf[:c], lineSep)

		switch {
		case err == io.EOF:
			return count

		case err != nil:
			return count
		}
	}
}
