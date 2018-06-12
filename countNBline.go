package main

import (
	"os"
	// "fmt"
)

var FILENAME string
var MAX int

func countLine(fname string) int {
	f_stat, err := os.Stat(fname)
	check(err)

	size := f_stat.Size()
	scanner := return_reader_for_bzipfile(fname)
	writer := return_writer_for_bzipfile("tmp.bz2")

	for i:=0;i<2000;i++ {
		scanner.Scan()
		writer.Write([]byte(scanner.Text()))
	}
	writer.Close()

	f_stat_writer, err := os.Stat("tmp.bz2")
	check(err)

	size_writer := f_stat_writer.Size()

	// fmt.Printf("\nsize %d\n", size)
	// fmt.Printf("\nsize_writer %d\n", size_writer)
	exe_cmd("rm tmp.bz2")

	return (int(size) / int(size_writer)) * int(2000)

}
