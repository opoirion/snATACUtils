package main


import (
	"path"
	"fmt"
	"strings"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"bufio"
	"os"
	"log"
	"sync"
	"bytes"
	"io"
	"time"
)

const BUFFERSIZE = 100000 - 100000 % 4

var BUFFERR1 [BUFFERSIZE]string
var BUFFERR2 [BUFFERSIZE]string
var BUFFERI1 [BUFFERSIZE]string
var COUNTS [3]int

var MUTEX sync.Mutex


/*FormatingR1R2FastqUsingI1Only format R1 and R2 fastq files using a 10x I1 fastq index file */
func FormatingR1R2FastqUsingI1Only(filenameR1 string, filenameR2 string, filenameI1 string) {
	MUTEX = sync.Mutex{}
	tStart := time.Now()

	var scannerI1 * bufio.Scanner
	var scannerR1 * bufio.Scanner
	var scannerR2 * bufio.Scanner

	var fileI1 * os.File
	var fileR1 * os.File
	var fileR2 * os.File

	var waiting sync.WaitGroup
	var begining, end, i int

	if BUFFERSIZE % 4 != 0 {
		log.Fatal("ERROR %d needs to be a multiple of 4 !", BUFFERSIZE)
	}

	chunk := BUFFERSIZE / 4 / NB_THREADS
	chunk = chunk - chunk % 4

	ext := path.Ext(filenameR1)
	outFilenameR1 := fmt.Sprintf("%s%s%s.demultiplexed.R1.fastq%s",
		OUTPUT_PATH,
		strings.TrimSuffix(filenameR1, fmt.Sprintf(".fastq%s", ext)),
		OUTPUT_TAG_NAME, ext)

	outFilenameR2 := fmt.Sprintf("%s%s%s.demultiplexed.R1.fastq%s",
		OUTPUT_PATH,
		strings.TrimSuffix(filenameR2, fmt.Sprintf(".fastq%s", ext)),
		OUTPUT_TAG_NAME, ext)

	writerR1 := utils.ReturnWriter(outFilenameR1)
	writerR2 := utils.ReturnWriter(outFilenameR2)

	defer writerR1.Close()
	defer writerR2.Close()

	scannerR1, fileR1 = utils.ReturnReader(filenameR1, 0, false)
	scannerR2, fileR2 = utils.ReturnReader(filenameR2, 0, false)
	scannerI1, fileI1 = utils.ReturnReader(filenameI1, 0, false)

	defer fileR1.Close()
	defer fileR2.Close()
	defer fileI1.Close()


	mainloop:
	for {
		waiting.Add(3)

		go fillBuffer(scannerR1, &BUFFERR1, &waiting, 0)
		go fillBuffer(scannerR2, &BUFFERR2, &waiting, 1)
		go fillBuffer(scannerI1, &BUFFERI1, &waiting, 2)

		waiting.Wait()

		switch {
		case COUNTS[0] != COUNTS[1] || COUNTS[0] != COUNTS[2] || COUNTS[1] != COUNTS[2]:
			log.Fatal(fmt.Sprintf("!!!! Error: Different read counts: %v found for the input fastq files!",
				COUNTS))
		case COUNTS[0] % 4 != 0:
			log.Fatal(fmt.Sprintf("!!!! Error: Number of lines not a multiple of 4: %d\n", COUNTS[0]))

		case COUNTS[0] == 0:
			break mainloop
		}

		end = chunk
		begining = 0

		for i = 0; i < NB_THREADS-1; i++ {
			if end >= COUNTS[0] {
				break
			}

			waiting.Add(1)
			go writeOutputFastq(begining, end, &writerR1, &writerR2, &waiting)

			end += chunk
			begining += chunk
		}

		waiting.Add(1)
		go writeOutputFastq(begining, COUNTS[0], &writerR1, &writerR2, &waiting)

		waiting.Wait()
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Formatting R1 R2 fastq files done in time: %f s \n", tDiff.Seconds())
	fmt.Printf("file: %s created!\n file: %s created!\n", outFilenameR1, outFilenameR2)

}

func fillBuffer(scanner * bufio.Scanner, buffer * [BUFFERSIZE]string, waiting * sync.WaitGroup, countPos int) {
	defer waiting.Done()
	count := 0

	for scanner.Scan() {
		buffer[count] = scanner.Text()
		count++

		if count >= BUFFERSIZE {
			break
		}
	}

	COUNTS[countPos] = count
}

func writeOutputFastq(begining, end int, writerR1, writerR2 *io.WriteCloser, waiting * sync.WaitGroup) {
	defer waiting.Done()

	var bufferR1 bytes.Buffer
	var bufferR2 bytes.Buffer

	defer bufferR1.Reset()
	defer bufferR2.Reset()

	if (end - begining) % 4 != 0 {
		log.Fatal(
			fmt.Sprintf("!!!! ERROR Number of lines  end: %d and begining: %d  needs to be modulo 3!", end, begining))
	}

	count := 4

	for i := begining; i < end;i++ {
		switch{
		case count == 4:
			if BUFFERR1[i][0] != '@' || BUFFERR2[i][0] != '@' || BUFFERI1[i][0] != '@' {
				log.Fatal(fmt.Sprintf("Error reads R1 %s and R2 %s I1 %s not sync for i %d\n",
					BUFFERR1[i], BUFFERR2[i], BUFFERI1[i], i))
			}

			bufferR1.WriteRune('@')
			bufferR1.WriteString(BUFFERI1[i+1])
			bufferR1.WriteRune(':')
			bufferR1.WriteString(BUFFERR1[i][1:])
			bufferR1.WriteRune('\n')

			bufferR2.WriteRune('@')
			bufferR2.WriteString(BUFFERI1[i+1])
			bufferR2.WriteRune(':')
			bufferR2.WriteString(BUFFERR2[i][1:])
			bufferR2.WriteRune('\n')

			count = 1

		default:
			bufferR1.WriteString(BUFFERR1[i])
			bufferR1.WriteRune('\n')
			bufferR2.WriteString(BUFFERR2[i])
			bufferR2.WriteRune('\n')

			count++
		}
	}

	MUTEX.Lock()
	(*writerR1).Write(bufferR1.Bytes())
	(*writerR2).Write(bufferR2.Bytes())
	MUTEX.Unlock()
}
