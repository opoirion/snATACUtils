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

//SymbolType Descriptor of a symbol: either a string or int slice
type SymbolType struct {
	SymbolPos []int
	SymbolStr string
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

//StringToPeakWithPos Convert Peak string to peak
func (peak * Peak) StringToPeakWithPos(str string, refPos [3]int) {
	var err1, err2 error
	split := strings.Split(str, "\t")

	(*peak).Start, err1 = strconv.Atoi(split[refPos[1]])
	(*peak).End, err2 = strconv.Atoi(split[refPos[2]])

	if err1 != nil || err2 != nil {
		panic(fmt.Sprintf(
			"Error when converting Peak: %s cannot be used as int ####\n",
			str))
	}

	(*peak).Slice[0] = split[refPos[0]]
	(*peak).Slice[1] = split[refPos[1]]
	(*peak).Slice[2] = split[refPos[2]]
}

//StringToPeakWithPosAndStart Convert Peak string to peak
func (peak * Peak) StringToPeakWithPosAndStart(str string, refPosList []int, start int) {
	var refPos [3]int

	if len(refPosList) < start + 3 {
		panic(fmt.Sprintf("Size error with refPosList: %d and start %d",
			refPosList, start))
	}

	refPos[0] = refPosList[0 + start]
	refPos[1] = refPosList[1 + start]
	refPos[2] = refPosList[2 + start]

	(*peak).StringToPeakWithPos(str, refPos)
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
var PEAKSYMBOLDICT map[Peak][]string


/*LoadSymbolFile  peaksymbolfile, peakfile  Filename*/
func LoadSymbolFile(peaksymbolfile, peakfile  Filename) {
	var scannerPeak *bufio.Scanner
	var filePeak *os.File
	var split []string
	var peakl Peak
	var symbol string

	PEAKSYMBOLDICT = make(map[Peak][]string)

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

		PEAKSYMBOLDICT[peakl] = append(PEAKSYMBOLDICT[peakl], symbol)
	}
}

/*LoadRefBedFileWithSymbol  peaksymbolfile, peakfile  Filename*/
func LoadRefBedFileWithSymbol(peaksymbolfile Filename) {
	symbol := SymbolType{}
	symbol.SymbolPos = []int{3}
	loadRefBedFileWithSymbol(peaksymbolfile, "\t", symbol, []int{0, 1, 2})
}

/*LoadRefCustomFileWithSymbol  peaksymbolfile, peakfile  Filename*/
func LoadRefCustomFileWithSymbol(
	peaksymbolfile Filename,
	sep string,
	symbol SymbolType,
	refPos []int) {
	loadRefBedFileWithSymbol(peaksymbolfile, sep, symbol, refPos)
}

/*CheckIfPeakPosIsMutltipleOf3 check if list is multiple of 3 */
func CheckIfPeakPosIsMutltipleOf3(peakPos []int) (numberOfPeaks int) {
	if len(peakPos) % 3 != 0 {
		panic(fmt.Sprintf(
			"peakPos %d from should b a multiple of 3",
			peakPos))
	}

	numberOfPeaks = len(peakPos) / 3

	return numberOfPeaks
}

/*loadRefBedFileWithSymbol  peaksymbolfile, peakfile  Filename*/
func loadRefBedFileWithSymbol(
	peaksymbolfile Filename, sep string, symbol SymbolType, peakPos []int) {
	var peakl Peak
	var symbolSlice, split, peaksplit []string
	var pos, i int
	var symbolStr string
	var peakPosTriplet [3]int

	symbolSlice = make([]string, len(symbol.SymbolPos))
	peaksplit = make([]string, 3)
	PEAKSYMBOLDICT = make(map[Peak][]string)

	maxPeakPos := MaxIntList(append(peakPos[:], symbol.SymbolPos...))

	numberOfPeaks := CheckIfPeakPosIsMutltipleOf3(peakPos)

	scanner, file := peaksymbolfile.ReturnReader(0)
	defer CloseFile(file)

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), sep)

		if split[0][0] == '#' {
			continue
		}

		if len(split) < maxPeakPos {
			panic(fmt.Sprintf(
				"Error line %s from symbol file should have at least enough fields as decribe in pos index: %d\n",
				split, peakPos))
		}

		for peakNb := 0; peakNb < numberOfPeaks; peakNb++ {
			peakPosTriplet[0] = peakPos[0 + 3 * peakNb]
			peakPosTriplet[1] = peakPos[1 + 3 * peakNb]
			peakPosTriplet[2] = peakPos[2 + 3 * peakNb]

			for i, pos = range symbol.SymbolPos {
				symbolSlice[i] = split[pos]
			}

			for i, pos = range peakPosTriplet {
				peaksplit[i] = split[pos]
			}

			peakl.SplitToPeak(peaksplit)

			if symbol.SymbolStr != "" {
				symbolStr = symbol.SymbolStr
			} else {
				symbolStr = strings.Join(symbolSlice, sep)
			}

			PEAKSYMBOLDICT[peakl] = append(PEAKSYMBOLDICT[peakl], symbolStr)
		}
	}
}

/*CreatePeakIntervalTreeCustom ...*/
func CreatePeakIntervalTreeCustom(peakPos []int, sep string) {
	createPeakIntervalTree(peakPos, sep, false)
}


/*CreatePeakIntervalTree ...*/
func CreatePeakIntervalTree() {
	createPeakIntervalTree([]int{0, 1, 2}, "\t", false)
}

/*createPeakIntervalTree ...*/
func createPeakIntervalTree(peakPos []int, sep string, verbose bool) {
	var split []string
	var chroStr string

	var start, end int
	var err error
	var isInside bool

	tStart := time.Now()

	CHRINTERVALDICT = make(map[string]*interval.IntTree)
	INTERVALMAPPING = make(map[uintptr]string)

	numberOfPeaks := CheckIfPeakPosIsMutltipleOf3(peakPos)
	maxPeakPos := MaxIntList(peakPos)

	for key, pos := range PEAKIDDICT {
		split = strings.Split(key, sep)

		if len(split) < maxPeakPos {
			panic(fmt.Sprintf(
				"Error from createPeakIntervalTree. Line %s from symbol file should have at least enough fields as decribe in pos index: %d\n",
				split, peakPos))
		}

		for peakNb := 0; peakNb < numberOfPeaks; peakNb++ {

			chroStr = split[peakPos[0 + 3 * peakNb]]

			start, err = strconv.Atoi(split[peakPos[1 + 3 * peakNb]])
			Check(err)

			end, err = strconv.Atoi(strings.Trim(split[peakPos[2 + 3 * peakNb]], "\n"))
			Check(err)

			inter := IntInterval{
				Start: start, End: end}
			inter.UID = uintptr(uintptr(pos))

			if _, isInside = CHRINTERVALDICT[chroStr];!isInside {
				CHRINTERVALDICT[chroStr] = &interval.IntTree{}
			}

			err = CHRINTERVALDICT[chroStr].Insert(inter, false)
			Check(err)

			INTERVALMAPPING[inter.ID()] = key
		}
	}

	tDiff := time.Since(tStart)

	if verbose {
		fmt.Printf("Create peak index done in time: %f s \n", tDiff.Seconds())
	}
}


/*createPeakIntervalTreeObject create a peak intervall dict object*/
func createPeakIntervalTreeObject(peakiddict map[string]uint, peakPos []int, verbose bool) (
	intervalObject PeakIntervalTreeObject) {

	var chroStr string
	var err error
	var isInside bool
	var peak Peak
	var peakPosTriplet [3]int

	numberOfPeaks := CheckIfPeakPosIsMutltipleOf3(peakPos)

	tStart := time.Now()

	intervalObject.Chrintervaldict = make(map[string]*interval.IntTree)
	intervalObject.Intervalmapping = make(map[uintptr]string)
	intervalObject.Peakiddict = &peakiddict

	for key, pos := range peakiddict {
		for peakNb := 0; peakNb < numberOfPeaks; peakNb++ {
			peakPosTriplet[0] = peakPos[0 + 3 * peakNb]
			peakPosTriplet[1] = peakPos[1 + 3 * peakNb]
			peakPosTriplet[2] = peakPos[2 + 3 * peakNb]

			peak.StringToPeakWithPos(key, peakPosTriplet)
			chroStr = peak.Chr()

			int := IntInterval{
				Start: peak.Start, End: peak.End}
			int.UID = uintptr(uintptr(pos))

			if _, isInside = intervalObject.Chrintervaldict[chroStr];!isInside {
				intervalObject.Chrintervaldict[chroStr] = &interval.IntTree{}
			}

			err = intervalObject.Chrintervaldict[chroStr].Insert(int, false)
			Check(err)

			intervalObject.Intervalmapping[int.ID()] = peak.PeakToString()
		}
	}

	tDiff := time.Since(tStart)

	if verbose {
		fmt.Printf("Create peak index done in time: %f s \n", tDiff.Seconds())
	}

	return intervalObject
}


/*CreatePeakIntervalTreeObjectFromFile create a peak intervall dict object*/
func CreatePeakIntervalTreeObjectFromFile(bedfile Filename, sep string, peakPos []int) (
	intervalObject PeakIntervalTreeObject) {

	peakiddict := LoadPeaksDictCustom(bedfile, sep, peakPos)

	intervalObject = createPeakIntervalTreeObject(peakiddict, peakPos, false)
	return intervalObject
}


/*LoadPeaksDict load peak file return map[string]int*/
func LoadPeaksDict(fname Filename) (peakiddict map[string]uint)  {
	peakiddict = make(map[string]uint)

	loadPeaks(fname, peakiddict, "\t", []int{0, 1, 2})

	return peakiddict
}


/*LoadPeaksDictCustom load peak file return map[string]int*/
func LoadPeaksDictCustom(fname Filename, sep string, peakPos []int) (
	peakiddict map[string]uint)  {
	peakiddict = make(map[string]uint)

	loadPeaks(fname, peakiddict, sep, peakPos)

	return peakiddict
}

/*LoadPeaks load peak file globally*/
func LoadPeaks(fname Filename) int {
	PEAKIDDICT = make(map[string]uint)

	return loadPeaks(fname, PEAKIDDICT, "\t", []int{0, 1, 2})
}

/*LoadPeaksCustom load peak file globally*/
func LoadPeaksCustom(fname Filename, sep string, peakPos []int) int {
	PEAKIDDICT = make(map[string]uint)

	return loadPeaks(fname, PEAKIDDICT, sep, peakPos)
}

/*loadPeaks load peak file globally*/
func loadPeaks(fname Filename, peakiddict map[string]uint, sep string, peakPos []int) int {
	var scanner *bufio.Scanner
	var file *os.File
	var line string

	scanner, file = fname.ReturnReader(0)

	defer CloseFile(file)

	count := uint(0)
	max := MaxIntList(peakPos)
	nbPeaks := CheckIfPeakPosIsMutltipleOf3(peakPos)

	for scanner.Scan() {
		line = scanner.Text()


		if line[0] == '#' {
			continue
		}

		checkIfLineCanBeSplitIntoPeaks(line, sep, peakPos, max, nbPeaks)

		peakiddict[line] = count
		count++
	}

	return int(count)
}


func checkIfLineCanBeSplitIntoPeaks(line, sep string, peakPos []int, peakMax, nbPeaks int) {
	var err1, err2 error

	split := strings.Split(line, sep)

	if len(split) < peakMax {
		panic(fmt.Sprintf(
			"line: %s cannot be splitted in more than %d with separator: %s to match peak position: %d",
			line, peakMax, sep, peakPos))
	}

	for peakNb := 0; peakNb < nbPeaks; peakNb++ {
		_, err1 = strconv.Atoi(split[1 + 3 * peakNb])
		_, err2 = strconv.Atoi(split[2 + 3 * peakNb])

		if err1 != nil || err2 != nil {
			panic(fmt.Sprintf(
				"line: %s cannot be converted into peak Position: %d",
				line, peakPos))
		}
	}

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


/*MaxIntList int give the max of list */
func MaxIntList(intlist []int) (max int) {

	for _, pos := range intlist {
		if pos > max {
			max = pos
		}
	}

	return max
}
