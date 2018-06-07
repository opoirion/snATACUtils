package main

import "fmt"
import "flag"
import "os"
import "log"
import "bufio"
import "strings"


var filename string
var taglength int

func main() {
	counter := make(map[string]int)
	flag.StringVar(&filename, "fastq", "", "fastq file name")
	flag.Parse()

	_, err := os.Stat(filename)

	fastq, err := os.Open(filename)

	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("fastq file analyzed: %s\n", filename)

	scanner := bufio.NewScanner(fastq)

	for scanner.Scan() {
		line := scanner.Text()

		if line[0] != '@' {
			continue
		}

		tag := strings.Split(line, ":")[0][1:]

		if taglength == 0 {
			taglength = len(tag)
		}

		if len(tag) != taglength {
			continue
		}

		fmt.Printf("%s\n", tag)
		counter[tag] += 1

	}

	fmt.Printf("length: %d\n", len(counter))

	for key, value := range counter {
		fmt.Printf("%s %d %d\n", key, value, len(key))

	}

}
