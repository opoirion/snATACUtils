package atacdemultiplexutils


import (
	"bufio"
	"os"
	"strings"
	"github.com/biogo/store/interval"
	"time"
	"fmt"
	"strconv"
	"github.com/jinzhu/copier"
)


//IntInterval Integer-specific intervals
type IntInterval struct {
	Start, End int
	UID        uintptr
	Payload    interface{}
}


//Overlap rule for two Interval
func (i IntInterval) Overlap(b interval.IntRange) bool {
	// Search for intersection
	return i.End >= b.Start && i.Start <= b.End
}

//ID Return the ID of Interval
func (i IntInterval) ID() uintptr {
	return i.UID
}
//Range Return the range of Interval
func (i IntInterval) Range() interval.IntRange {
	return interval.IntRange{i.Start, i.End}
}

//String Return the string re[ of Interval
func (i IntInterval) String() string {
	return fmt.Sprintf("(%d, %d) id: %d ####\n", i.Start, i.End, i.ID())
}


/*PEAKIDDICT peak ID<->pos */
var PEAKIDDICT map[string]uint

/*CHRINTERVALDICT chr ID <-> interval tree */
var CHRINTERVALDICT map[string]*interval.IntTree

/*CHRINTERVALDICTTHREAD threadNB -> chr ID -> pos */
var CHRINTERVALDICTTHREAD map[int]map[string]*interval.IntTree

/*INTERVALMAPPING peak ID pos <->pos */
var INTERVALMAPPING map[uintptr]string


/*CreatePeakIntervalTree ...*/
func CreatePeakIntervalTree() {
	var split []string
	var chroStr string

	var start, end int
	var err error
	var isInside bool

	fmt.Printf("create peak interval tree...\n")
	tStart := time.Now()

	CHRINTERVALDICT = make(map[string]*interval.IntTree)
	INTERVALMAPPING = make(map[uintptr]string)

	for key, pos := range PEAKIDDICT {
		split = strings.Split(key, "\t")
		chroStr = split[0]

		start, err = strconv.Atoi(split[1])
		Check(err)

		end, err = strconv.Atoi(strings.Trim(split[2], "\n"))
		Check(err)

		int := IntInterval{
			Start: start, End: end}
		int.UID = uintptr(uintptr(pos))

		if _, isInside = CHRINTERVALDICT[chroStr];!isInside {
			CHRINTERVALDICT[chroStr] = &interval.IntTree{}
		}

		err = CHRINTERVALDICT[chroStr].Insert(int, false)
		Check(err)

		INTERVALMAPPING[int.ID()] = key
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Create peak index done in time: %f s \n", tDiff.Seconds())
}



/*LoadPeaks load peak file and return peak peak id -> dict*/
func LoadPeaks(fname Filename) int {
	var scanner *bufio.Scanner
	var file *os.File
	var split []string
	var err1, err2 error
	var line string

	scanner, file = fname.ReturnReader(0)

	defer CloseFile(file)


	var count uint

	PEAKIDDICT = make(map[string]uint)
	count = 0

	for scanner.Scan() {
		line = scanner.Text()

		split = strings.Split(line, "\t")

		if len(split) < 3 {
			panic(fmt.Sprintf(
				"Peak: %s at line %d from file %s cannot be cut in chr int int ####\n",
				line, count, fname))
		}

		_, err1 = strconv.Atoi(split[1])
		_, err2 = strconv.Atoi(split[2])

		if err1 != nil || err2 != nil {
			panic(fmt.Sprintf(
				"Peak positions : %s at line %d from file %s cannot be used as int ####\n",
				line, count, fname))
		}

		PEAKIDDICT[line] = count
		count++
	}

	return int(count)
}


/*LoadPeaksAndTrim load peak file and return peak peak id trimmed for "chr" -> dict*/
func LoadPeaksAndTrim(fname Filename) int {
	var scanner *bufio.Scanner
	var file *os.File
	var line string
	var split []string
	var err1, err2 error

	scanner, file = fname.ReturnReader(0)

	defer CloseFile(file)


	var count uint

	PEAKIDDICT = make(map[string]uint)
	count = 0

	for scanner.Scan() {
		line = strings.TrimPrefix(scanner.Text(), "chr")
		split = strings.Split(line, "\t")

		if len(split) < 3 {
			panic(fmt.Sprintf(
				"Peak: %s at line %d from file %s cannot be cut in chr int int ####\n",
				line, count, fname))
		}

		_, err1 = strconv.Atoi(split[1])
		_, err2 = strconv.Atoi(split[2])

		if err1 != nil || err2 != nil {
			panic(fmt.Sprintf(
				"Peak positions : %s at line %d from file %s cannot be used as int ####\n",
				line, count, fname))
		}

		PEAKIDDICT[line] = count
		count++
	}

	return int(count)
}


/*LoadPeaksSubset load peak file  but using only a subset of peaks and return peak peak id -> dict*/
func LoadPeaksSubset(fname Filename, firstPeak, lastPeak int) {
	var scanner *bufio.Scanner
	var file *os.File

	peaknb := -1

	scanner, file = fname.ReturnReader(0)

	defer CloseFile(file)


	var count uint

	PEAKIDDICT = make(map[string]uint)
	count = 0

	for scanner.Scan() {
		peaknb++

		if peaknb < firstPeak {
			continue
		}

		if peaknb >= lastPeak {
			break
		}

		line := scanner.Text()
		PEAKIDDICT[line] = count
		count++
	}
}


/*InitIntervalDictsThreading Init interval dict threading map by copying the interval map for each trheads*/
func InitIntervalDictsThreading(threadnb int) {
       CHRINTERVALDICTTHREAD = make(map[int]map[string]*interval.IntTree)

       for i := 0;i< threadnb;i++ {
               CHRINTERVALDICTTHREAD[i] = make(map[string]*interval.IntTree)

               for key, tree := range CHRINTERVALDICT {
                       CHRINTERVALDICTTHREAD[i][key] = &interval.IntTree{}
                       err := copier.Copy(CHRINTERVALDICTTHREAD[i][key], tree)
		       Check(err)
               }
       }
}
