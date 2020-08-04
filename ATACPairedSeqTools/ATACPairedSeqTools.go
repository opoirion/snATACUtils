package main


import (
	"fmt"
	"flag"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"os"
	"path"
	"bytes"
	"sync"
	"io"
	"strings"
	"time"
)


/*FASTQPATH input fastq file*/
var FASTQPATH utils.Filename

/*OUTFILE output file*/
var OUTFILE string

/*FORMATNAME format name action*/
var FORMATNAME bool

/*SAMTOFASTQ sam to fastq*/
var SAMTOFASTQ bool

/*SAMFILE sam input*/
var SAMFILE utils.Filename

/*WAITING waiting group*/
var WAITING sync.WaitGroup

/*MUTEX1 mutex */
var MUTEX1 sync.Mutex

/*MUTEX2 mutex */
var MUTEX2 sync.Mutex

/*MUTEXW mutex */
var MUTEXW sync.Mutex

/*REPLACEINPUT  edit input fastq file instrad of creating a new file */
var REPLACEINPUT bool


func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `
#### ATACPairesdSeqTools USAGE

ATACPairesdSeqTools -format -in <fastq file> -out <fastq file>
`)
		flag.PrintDefaults()
	}

	flag.StringVar(&OUTFILE, "out", "", "output fastq file ")
	flag.Var(&FASTQPATH, "fastq", "input fastq file ")
	flag.Var(&SAMFILE, "sam", "input sam file ")
	flag.BoolVar(&FORMATNAME, "format", false, "Format fastq file by integrating name")
	flag.BoolVar(&SAMTOFASTQ, "sam_to_fastq", false, "Format SAM to FASTQ")
	flag.BoolVar(&REPLACEINPUT, "edit", false,
		`edit input fastq file instead of creating a new file`)

	flag.Parse()

	tStart := time.Now()
	switch {
	case FORMATNAME:
		formatName()
	case SAMTOFASTQ:
		samToFastq()
	default:
		flag.Usage()
	}

	tDiff := time.Since(tStart)
	fmt.Printf("done in time: %f s \n", tDiff.Seconds())

	if REPLACEINPUT {
		utils.Check(os.Remove(string(FASTQPATH)))
		utils.Check(os.Rename(OUTFILE, string(FASTQPATH)))
		fmt.Printf("File: %s edited\n", FASTQPATH)
	} else {
		fmt.Printf("File: %s created\n", OUTFILE)
	}
}


func samToFastq() {
	if SAMFILE == "" {
		fmt.Printf("#### Error -sam should be set!\n")
		flag.Usage()
		os.Exit(1)
	}

	if OUTFILE == "" {
		ext := path.Ext(SAMFILE.String())

		OUTFILE = fmt.Sprintf("%s.formatted.fastq.gz",
			SAMFILE[:len(SAMFILE) - len(ext)])
	}

	var split, split2 []string
	var line string
	var splitl int

	samReader, file := utils.ReturnReader(SAMFILE.String(), 0)
	defer utils.CloseFile(file)

	writer := utils.ReturnWriter(OUTFILE)
	defer utils.CloseFile(writer)

	buffer1, buffer2  := bytes.Buffer{}, bytes.Buffer{}
	MUTEX1, MUTEX2, MUTEXW = sync.Mutex{}, sync.Mutex{}, sync.Mutex{}
	WAITING = sync.WaitGroup{}

	isBuffer1 := true
	mutex := &MUTEX1
	buffer := &buffer1

	count := 0
	mutex.Lock()

	for samReader.Scan() {
		line = samReader.Text()

		if line[0] == '@' {
			continue
		}

		split2 = strings.Split(line, "\t")

		split = strings.Split(split2[0], ":")
		splitl = len(split)

		if split2[2] == "*" {
			continue
		}

		buffer.WriteRune('@')
		buffer.WriteString(split[0])
		buffer.WriteString(split2[9])
		buffer.WriteRune(':')
		buffer.WriteString(strings.Join(split[1:splitl-2], ":"))
		buffer.WriteRune(':')
		buffer.WriteString(split[0])
		buffer.WriteRune(':')
		buffer.WriteString(split2[2])
		buffer.WriteRune('\n')
		buffer.WriteString(split[splitl-2])
		buffer.WriteString("\n+\n")
		buffer.WriteString(split[splitl-1])
		buffer.WriteRune('\n')

		count++

		if count > 10000 {
			WAITING.Add(1)
			mutex.Unlock()
			go writeBuffer(buffer, &writer, mutex)

			if isBuffer1 {
				mutex = &MUTEX2
				buffer = &buffer2
				isBuffer1 = false
			} else {
				mutex = &MUTEX1
				buffer = &buffer1
				isBuffer1 = true
			}

			mutex.Lock()
		}
	}

	WAITING.Add(1)
	mutex.Unlock()
	writeBuffer(buffer, &writer, mutex)

	WAITING.Wait()
}


func formatName() {
	if FASTQPATH == "" {
		fmt.Printf("#### Error -fastq should be set!\n")
		flag.Usage()
		os.Exit(1)
	}

	var line string
	var splitl int
	var split[]string

	if OUTFILE == "" {
		ext := path.Ext(FASTQPATH.String())

		OUTFILE = fmt.Sprintf("%s.formatted%s",
			FASTQPATH[:len(FASTQPATH) - len(ext)], ext)
	}

	scanner, file := FASTQPATH.ReturnReader(0)
	defer utils.CloseFile(file)

	writer := utils.ReturnWriter(OUTFILE)
	defer utils.CloseFile(writer)

	buffer1, buffer2  := bytes.Buffer{}, bytes.Buffer{}
	MUTEX1, MUTEX2, MUTEXW = sync.Mutex{}, sync.Mutex{}, sync.Mutex{}
	WAITING = sync.WaitGroup{}

	isBuffer1 := true
	mutex := &MUTEX1
	buffer := &buffer1

	count := 0
	mutex.Lock()

	for scanner.Scan() {
		line = scanner.Text()

		if line[0] != '@' {
			buffer.WriteString(line)
			buffer.WriteRune('\n')
			count++
		} else {
			split = strings.Split(line, ":")
			splitl = len(split)
			buffer.WriteString(split[0])
			buffer.WriteRune('-')
			buffer.WriteString(split[splitl-4])
			buffer.WriteRune('-')
			buffer.WriteString(split[splitl-3])
			buffer.WriteRune('-')
			buffer.WriteString(split[splitl-2])
			buffer.WriteRune(':')
			buffer.WriteString(strings.Join(split[1:], ":"))
			buffer.WriteRune('\n')
			count++
		}

		if count > 1000000 {
			WAITING.Add(1)
			mutex.Unlock()
			go writeBuffer(buffer, &writer, mutex)

			if isBuffer1 {
				mutex = &MUTEX2
				buffer = &buffer2
				isBuffer1 = false
			} else {
				mutex = &MUTEX1
				buffer = &buffer1
				isBuffer1 = true
			}

			mutex.Lock()
		}
	}

	WAITING.Add(1)
	mutex.Unlock()
	writeBuffer(buffer, &writer, mutex)

	WAITING.Wait()
}

func writeBuffer(buffer * bytes.Buffer,
	writer *io.WriteCloser,
	mutex * sync.Mutex) {
	mutex.Lock()
	defer WAITING.Done()
	MUTEXW.Lock()
	_, err := (*writer).Write(buffer.Bytes())
	MUTEXW.Unlock()
	utils.Check(err)
	buffer.Reset()
	mutex.Unlock()
}
