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


//PeakIntervalTreeObject Peak IntervalTree Object
type PeakIntervalTreeObject struct {
	Chrintervaldict map[string]*interval.IntTree
	Intervalmapping map[uintptr]string
	Peakiddict *map[string]uint
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

//Peak Descriptor of a peak as string slice
type Peak struct{
	Slice [3]string
	Start, End int
}

//StringToPeak Convert Peak string to peak
func (peak * Peak) StringToPeak(str string) {
	var err1, err2 error
	split := strings.Split(str, "\t")

	(*peak).Start, err1 = strconv.Atoi(split[1])
	(*peak).End, err2 = strconv.Atoi(split[2])

	if err1 != nil || err2 != nil {
		panic(fmt.Sprintf(
			"Error when converting Peak: %s cannot be used as int ####\n",
			str))
	}

	(*peak).Slice[0] = split[0]
	(*peak).Slice[1] = split[1]
	(*peak).Slice[2] = split[2]
}

//SplitToPeak Convert string split to peak
func (peak * Peak) SplitToPeak(split []string) {
	var err1, err2 error

	(*peak).Start, err1 = strconv.Atoi(split[1])
	(*peak).End, err2 = strconv.Atoi(split[2])

	if err1 != nil || err2 != nil {
		panic(fmt.Sprintf(
			"Error when converting Peak: %s cannot be used as int ####\n",
			split))
	}

	(*peak).Slice[0] = split[0]
	(*peak).Slice[1] = split[1]
	(*peak).Slice[2] = split[2]
}

/*PeakToString Convert Peak to string*/
func (peak * Peak) PeakToString() (peakstr string)  {
	return fmt.Sprintf("%s\t%s\t%s", (*peak).Slice[0],
		(*peak).Slice[1],
		(*peak).Slice[2])

}

/*Chr return the chromosome of the peak */
func (peak * Peak) Chr() (chr string)  {
	return (*peak).Slice[0]
}

//StringToPeakNoCheck Convert Peak string to peak
func (peak * Peak) StringToPeakNoCheck(str string) {
	split := strings.Split(str, "\t")

	(*peak).Slice[0] = split[0]
	(*peak).Slice[1] = split[1]
	(*peak).Slice[2] = split[2]
}

/*PEAKIDDICT peak ID<->pos */
var PEAKIDDICT map[string]uint

/*CHRINTERVALDICT chr ID <-> interval tree */
var CHRINTERVALDICT map[string]*interval.IntTree

/*CHRINTERVALDICTTHREAD threadNB -> chr ID -> pos */
var CHRINTERVALDICTTHREAD map[int]map[string]*interval.IntTree

/*INTERVALMAPPING peak ID pos <->pos */
var INTERVALMAPPING map[uintptr]string

/*PEAKSYMBOLDICT map[peak]symbol */
var PEAKSYMBOLDICT map[Peak]string


/*LoadSymbolFile  peaksymbolfile, peakfile  Filename*/
func LoadSymbolFile(peaksymbolfile, peakfile  Filename) {
	var scannerPeak *bufio.Scanner
	var filePeak *os.File
	var split []string
	var peakl Peak
	var symbol string

	PEAKSYMBOLDICT = make(map[Peak]string)

	if peaksymbolfile == "" {
		return
	}

	isOption1 := true

	if peakfile == "" {
		isOption1 = false
	} else {
		scannerPeak, filePeak = peakfile.ReturnReader(0)
		defer CloseFile(filePeak)
	}

	scanner, file := peaksymbolfile.ReturnReader(0)
	defer CloseFile(file)

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), "\t")

		if len(split) == 4 {
			isOption1 = false
		}

		if !isOption1 && len(split) != 4 {
			panic(fmt.Sprintf(
				"Error line %s from symbol file should be <symbol>\t<chromosome>\t<start>\t<stop>\n",
				split))
		}

		symbol = split[0]

		if isOption1 {
			scannerPeak.Scan()
			peakl.StringToPeak(scannerPeak.Text())

		} else {
			peakl.SplitToPeak(split)
		}

		PEAKSYMBOLDICT[peakl] = symbol
	}
}

/*LoadRefBedFileWithSymbol  peaksymbolfile, peakfile  Filename*/
func LoadRefBedFileWithSymbol(peaksymbolfile Filename) {
	loadRefBedFileWithSymbol(peaksymbolfile, "\t", []int{3})
}

/*LoadRefCustomFileWithSymbol  peaksymbolfile, peakfile  Filename*/
func LoadRefCustomFileWithSymbol(peaksymbolfile Filename, sep string, symbolPos []int) {
	loadRefBedFileWithSymbol(peaksymbolfile, sep, symbolPos)
}

/*loadRefBedFileWithSymbol  peaksymbolfile, peakfile  Filename*/
func loadRefBedFileWithSymbol(peaksymbolfile Filename, sep string, symbolPos []int) {
	var split []string
	var peakl Peak
	var symbol []string
	var pos, i int

	symbol = make([]string, len(symbolPos))
	PEAKSYMBOLDICT = make(map[Peak]string)


	scanner, file := peaksymbolfile.ReturnReader(0)
	defer CloseFile(file)

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), sep)

		if len(split) < 4 {
			panic(fmt.Sprintf(
				"Error line %s from symbol file should be <chromosome>\t<start>\t<stop>\tsymbol>\n",
				split))
		}

		for i, pos = range symbolPos {
			symbol[i] = split[pos]
		}

		peakl.SplitToPeak(split)

		PEAKSYMBOLDICT[peakl] = strings.Join(symbol, sep)
	}
}

/*CreatePeakIntervalTreeCustom ...*/
func CreatePeakIntervalTreeCustom(peakPos [3]int, sep string) {
	createPeakIntervalTree(peakPos, sep)
}


/*CreatePeakIntervalTree ...*/
func CreatePeakIntervalTree() {
	createPeakIntervalTree([3]int{0, 1, 2}, "\t")
}

/*createPeakIntervalTree ...*/
func createPeakIntervalTree(peakPos [3]int, sep string) {
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
		split = strings.Split(key, sep)
		chroStr = split[peakPos[0]]

		start, err = strconv.Atoi(split[peakPos[1]])
		Check(err)

		end, err = strconv.Atoi(strings.Trim(split[peakPos[2]], "\n"))
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


/*CreatePeakIntervalTreeObject create a peak intervall dict object*/
func CreatePeakIntervalTreeObject(peakiddict map[string]uint) (
	intervalObject PeakIntervalTreeObject) {

	var split []string
	var chroStr string

	var start, end int
	var err error
	var isInside bool

	fmt.Printf("create peak interval tree...\n")
	tStart := time.Now()

	intervalObject.Chrintervaldict = make(map[string]*interval.IntTree)
	intervalObject.Intervalmapping = make(map[uintptr]string)
	intervalObject.Peakiddict = &peakiddict

	for key, pos := range peakiddict {
		split = strings.Split(key, "\t")
		chroStr = split[0]

		start, err = strconv.Atoi(split[1])
		Check(err)

		end, err = strconv.Atoi(strings.Trim(split[2], "\n"))
		Check(err)

		int := IntInterval{
			Start: start, End: end}
		int.UID = uintptr(uintptr(pos))

		if _, isInside = intervalObject.Chrintervaldict[chroStr];!isInside {
			intervalObject.Chrintervaldict[chroStr] = &interval.IntTree{}
		}

		err = intervalObject.Chrintervaldict[chroStr].Insert(int, false)
		Check(err)

		intervalObject.Intervalmapping[int.ID()] = key
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Create peak index done in time: %f s \n", tDiff.Seconds())

	return intervalObject
}

/*CreatePeakIntervalTreeObjectFromFile create a peak intervall dict object*/
func CreatePeakIntervalTreeObjectFromFile(bedfile Filename) (intervalObject PeakIntervalTreeObject) {
	peakiddict := LoadPeaksDict(bedfile)

	intervalObject = CreatePeakIntervalTreeObject(peakiddict)
	return intervalObject
}


/*LoadPeaksDict load peak file return map[string]int*/
func LoadPeaksDict(fname Filename) (peakiddict map[string]uint)  {
	peakiddict = make(map[string]uint)

	loadPeaks(fname, peakiddict, "\t", [3]int{0, 1, 2})

	return peakiddict
}

/*LoadPeaks load peak file globally*/
func LoadPeaks(fname Filename) int {
	PEAKIDDICT = make(map[string]uint)

	return loadPeaks(fname, PEAKIDDICT, "\t", [3]int{0, 1, 2})
}

/*LoadPeaksCustom load peak file globally*/
func LoadPeaksCustom(fname Filename, sep string, peakPos [3]int) int {
	PEAKIDDICT = make(map[string]uint)

	return loadPeaks(fname, PEAKIDDICT, sep, peakPos)
}

/*loadPeaks load peak file globally*/
func loadPeaks(fname Filename, peakiddict map[string]uint, sep string, pos [3]int) int {
	var scanner *bufio.Scanner
	var file *os.File
	var split []string
	var err1, err2 error
	var line string

	scanner, file = fname.ReturnReader(0)

	defer CloseFile(file)

	var count uint
	count = 0

	max := 3

	for p := range pos {
		if p >= max {
			max = p + 1
		}
	}

	for scanner.Scan() {
		line = scanner.Text()

		split = strings.Split(line, sep)

		if len(split) < max {
			panic(fmt.Sprintf(
				"Peak: %s at line %d from file %s cannot be cut in chr int int ####\n",
				line, count, fname))
		}

		_, err1 = strconv.Atoi(split[pos[1]])
		_, err2 = strconv.Atoi(split[pos[2]])

		if err1 != nil || err2 != nil {
			panic(fmt.Sprintf(
				"Peak positions : %s at line %d from file %s cannot be used as int ####\n",
				line, count, fname))
		}

		peakiddict[line] = count
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
