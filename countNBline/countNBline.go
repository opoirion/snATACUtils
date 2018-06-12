package main

import(
	"fmt";
	"flag";
	"bufio";
	"os";
	"log";
	"time";
)




var FILENAME string
var MAX int

func main() {
	flag.StringVar(&FILENAME, "filename", "", "name of the file to count the lines")
	flag.IntVar(&MAX, "max nb lines", 0, "max number of lines")
	flag.Parse()
	fmt.Printf("%s\n", FILENAME)
	t_start := time.Now()
	t_diff := time.Now().Sub(t_start)

	nbLines := countLine(FILENAME)

	fmt.Printf("number of lines: %d\n", nbLines)
	fmt.Printf("time: %f s\n", t_diff.Seconds())

}


func countLine(filename string) int {

	file, err := os.Open(filename)

	check(err)

	scanner := bufio.NewScanner(file)

	nbLines := 0

	for scanner.Scan() {
		nbLines++
	}

	return nbLines
}


func check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}
