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
	"strconv"
)

/*BUFFERSIZE buufer size used to parse strings from ungzipped file */
const BUFFERSIZE = 100000 - 100000 % 4

/*BUFFERR1 buffer 1 */
var BUFFERR1 [BUFFERSIZE]string

/*BUFFERR2 buffer 2 */
var BUFFERR2 [BUFFERSIZE]string

/*BUFFERI1 buffer index I1 */
var BUFFERI1 [BUFFERSIZE]string

/*COUNTS reqad count array for the 3 files */
var COUNTS [3]int

/*MUTEX main mutex */
var MUTEX sync.Mutex

/*CHANMUTEX Mutex to lock chan storing*/
var CHANMUTEX sync.Mutex

/*CHANWAITING waiting group to lock chan storing*/
var CHANWAITING sync.WaitGroup

/*LOGSDICTMAP global log dict file -> field type -> entity -> count  */
var LOGSDICTMAP map[string]map[string]map[string]int


/*FormatingR1R2FastqUsingI1Only format R1 and R2 fastq files using a 10x I1 fastq index file */
func FormatingR1R2FastqUsingI1Only(filenameR1 string, filenameR2 string, filenameI1 string) {
	MUTEX = sync.Mutex{}
	CHANMUTEX = sync.Mutex{}

	tStart := time.Now()

	LOGSDICTMAP = make(map[string]map[string]map[string]int)

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

	fmt.Printf("#### Launching ATACdemultiplex using only I1 #### \n")

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

	scannerR1, fileR1 = utils.ReturnReader(filenameR1, 0)
	scannerR2, fileR2 = utils.ReturnReader(filenameR2, 0)
	scannerI1, fileI1 = utils.ReturnReader(filenameI1, 0)

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

	CHANWAITING.Wait()

	go storeDictFromChan(LOG_CHAN, "stats")
	go storeDictFromChan(LOG_INDEX_CELL_CHAN, "barcodes")
	go storeDictFromChan(LOG_INDEX_READ_CHAN, "failed_barcodes")
	CHANWAITING.Wait()

	go writeReportFrom10x("stats")
	go writeReportFrom10x("barcodes")

	if !USENOINDEX {
		go writeReportFrom10x("failed_barcodes")
	}

	CHANWAITING.Wait()

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

	var successP7 bool
	// var replicateP7 int

	logsIndexCell := initLog(LOG_INDEX_TYPE)
	logsIndexCellFail := initLog(LOG_INDEX_TYPE)

	logs := initLog(LOG_TYPE)

	if (end - begining) % 4 != 0 {
		log.Fatal(
			fmt.Sprintf("!!!! ERROR Number of lines  end: %d and begining: %d  needs to be modulo 3!", end, begining))
	}

	count := 4
	success := false

loop:
	for i := begining; i < end;i++ {
		switch{
		case count == 4:
			if BUFFERR1[i][0] != '@' || BUFFERR2[i][0] != '@' || BUFFERI1[i][0] != '@' {
				log.Fatal(fmt.Sprintf("Error reads R1 %s and R2 %s I1 %s not sync for i %d\n",
					BUFFERR1[i], BUFFERR2[i], BUFFERI1[i], i))
			}

			if !USENOINDEX {
				successP7, BUFFERI1[i+1], _ = checkOneIndex(BUFFERI1[i+1], "p7")

				if !successP7 {
					if WRITE_LOGS {
						logsIndexCellFail[fmt.Sprintf("fail_p7")].dict[BUFFERI1[i+1]]++
						logs["stats"].dict["Number of reads"]++
						count++
						success = false
						continue loop
					}
				}
			}

			success = true

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

			if WRITE_LOGS {
				logsIndexCell[fmt.Sprintf("success_p7_repl1")].dict[BUFFERI1[i+1]]++
				logs["stats"].dict["Number of reads"]++
			}

		default:
			if success {
				bufferR1.WriteString(BUFFERR1[i])
				bufferR1.WriteRune('\n')
				bufferR2.WriteString(BUFFERR2[i])
				bufferR2.WriteRune('\n')
			}

			count++
		}
	}

	MUTEX.Lock()
	(*writerR1).Write(bufferR1.Bytes())
	(*writerR2).Write(bufferR2.Bytes())
	MUTEX.Unlock()

	go processLogChan(LOG_CHAN, logs, "stats")
	go processLogChan(LOG_INDEX_CELL_CHAN, logsIndexCell, "barcodes")

	if !USENOINDEX {
		go processLogChan(LOG_INDEX_READ_CHAN, logsIndexCellFail, "failed_barcodes")
	}
}


func processLogChan(channels  map[string]chan StatsDict,
	logs map[string]*StatsDict,
	logFileKey string) {
	CHANWAITING.Add(1)
	defer CHANWAITING.Done()
	for key, value := range logs {
		if len(channels[key]) == cap(channels[key]) {
			storeDictFromChan(channels, logFileKey)
		}

		channels[key] <- *value
	}

}

func storeDictFromChan(channels  map[string]chan StatsDict, logFileKey string) {

	var isInside bool
	CHANWAITING.Add(1)
	defer CHANWAITING.Done()

	CHANMUTEX.Lock()
	defer CHANMUTEX.Unlock()

	if _, isInside = LOGSDICTMAP[logFileKey];!isInside {
		LOGSDICTMAP[logFileKey] = make(map[string]map[string]int)
	}


	for typ, channel := range channels {
		if _, isInside = LOGSDICTMAP[logFileKey][typ];!isInside {
			LOGSDICTMAP[logFileKey][typ] = make(map[string]int)
		}

		loop:
		for {
			select{
			case statsDict := <-channel:
				dictit := statsDict.dict
				for key, value := range dictit {
					LOGSDICTMAP[logFileKey][typ][key] += value
				}
			default:
				break loop
			}
		}
	}
}


func writeReportFrom10x(logFileKey string) {
	if !WRITE_LOGS{
		return
	}
	tStart := time.Now()
	var buffer bytes.Buffer

	fmt.Printf("writing reports....\n")
	CHANWAITING.Add(1)
	defer CHANWAITING.Done()

	filename := fmt.Sprintf("%sreport%s_%s.log",
		OUTPUT_PATH, OUTPUT_TAG_NAME, logFileKey)

	file, err := os.Create(filename)
	Check(err)
	defer file.Close()

	for key, logs := range LOGSDICTMAP[logFileKey] {
		file.WriteString(fmt.Sprintf("#### %s\n", key))
		file.WriteString("#<key>\t<value>\n")

		switch SORT_LOGS {
		case true:
			rankedLogs := utils.RankByWordCountAndDeleteOldMap(&logs)

			for _, pair := range rankedLogs {
				buffer.WriteString(pair.Key)
				buffer.WriteRune('\t')
				buffer.WriteString(strconv.Itoa(pair.Value))
				buffer.WriteRune('\n')

				file.Write(buffer.Bytes())
				buffer.Reset()
			}

		default:
			for k, v := range logs {
				buffer.WriteString(k)
				buffer.WriteRune('\t')
				buffer.WriteString(strconv.Itoa(v))
				buffer.WriteRune('\n')

				file.Write(buffer.Bytes())
				buffer.Reset()
			}
		}
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Report written in: %f s \n", tDiff.Seconds())
}
