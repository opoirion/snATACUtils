package main

import (
	"flag"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"os"
	"fmt"
	"path"
	"github.com/biogo/store/interval"
	"strings"
	"strconv"
	"bytes"
	"time"
	"math"
	"io"
)



/*BEDFILENAME bed file name (input) */
var BEDFILENAME utils.Filename

/*REFBEDFILENAME bed file containing annotation as fourth column */
var REFBEDFILENAME utils.Filename

/*FILENAMEOUT  output file name output */
var FILENAMEOUT string

/*REPLACEINPUT  edit input bed file */
var REPLACEINPUT bool

/*IGNOREUNANNOATED ignore unnatotated peak */
var IGNOREUNANNOATED bool

/*WRITEINTERSECT write intersection only */
var WRITEINTERSECT bool

/*UNIQ write only unique output peaks */
var UNIQ bool

/*UNIQREF write only unique output peaks */
var UNIQREF bool

/*UNIQSYMBOL write only unique output peaks */
var UNIQSYMBOL bool

/*WRITEREF write_ref write bed region from reference file */
var WRITEREF bool

/*WRITEDIFF write only element that are not intersecting */
var WRITEDIFF bool

/*REFSEP separator used to identify the reference region in the -ref file */
var REFSEP string

/*REFPOS position of reference region in the -ref file */
var REFPOS string

/*SYMBOLPOS position of the coulmns used for annotations in the -ref file */
var SYMBOLPOS string

/*STDOUT write to stdout*/
var STDOUT bool


func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `
#################### MODULE TO ANNOTATE GENOMIC REGIONS FROM BED FILES ########################

"""Annotate bed file using a reference bed file containing the annotations"""
USAGE: ATACAnnotateRegions -bed <file> -ref <file> (optionnal -out <string> -unique -unique_ref -intersect -write_ref -edit -ref_sep "[3]int" -ref_symbol "[]int" -diff -stdout)

for -ref_sep and -ref_symbol options, the input should be a string of numbers separated by whitespace and delimited with ". -ref_sep needs exactly three positions: 1) for the chromosomes column, 2) for the begining and 3) for the end of the region

Example: ATACAnnotateRegions -bed regionToAnnotate.bed -ref referenceAnnotation.tsv -ref_sep "0 1 2" -ref_symbol "4 5"

Here the three first columns of referenceAnnotation.tsv will be used to identify chromosome (column 0), start (column 1), and end (column 2) of each region, and regionToAnnotate.bed will be annotatd using columns 4 and 5 from referenceAnnotation.tsv
`)
		flag.PrintDefaults()
	}

	flag.Var(&BEDFILENAME, "bed", "name of the bed file no annotate")
	flag.Var(&REFBEDFILENAME, "ref", "name of the reference bed file containing the annotations")
	flag.StringVar(&FILENAMEOUT, "out", "", "name the output file(s)")
	flag.BoolVar(&REPLACEINPUT, "edit", false, `edit input bed file instead of creating a new file`)
	flag.BoolVar(&IGNOREUNANNOATED, "ignore", false, `ignore unnatotated peak`)
	flag.BoolVar(&WRITEINTERSECT, "intersect", false, `write intersection only`)
	flag.BoolVar(&UNIQ, "unique", false, `write only unique output peaks`)
	flag.BoolVar(&UNIQREF, "unique_ref", false, `write only unique reference using the closest peak`)
	flag.BoolVar(&UNIQSYMBOL, "unique_symbols", true, `write only unique symbols per peak`)
	flag.BoolVar(&WRITEDIFF, "diff", false, `write bed region if no intersection is found`)
	flag.BoolVar(&WRITEREF, "write_ref", false, `write bed region from reference file`)
	flag.BoolVar(&STDOUT, "stdout", false, `write to stdout`)
	flag.StringVar(&REFSEP, "ref_sep", "\t", "separator to define the bed region for the ref file")
	flag.StringVar(&REFPOS, "ref_pos", "0 1 2", "separator to the bed region in ref for the ref file")
	flag.StringVar(&SYMBOLPOS, "symbol_pos", "3", "separator to the bed region in ref for the ref file")

	flag.Parse()

	tail := flag.Args()

	if len(tail) > 0 {
		panic(fmt.Sprintf(`Error wrongly formatted arguments: %s for -ref_sep and -ref_symbol options, the input should be a string of numbers separated by whitespace and delimited with ". -ref_sep needs exactly three positions: 1) for the chromosomes column, 2) for the begining and 3) for the end of the region

Example: ATACAnnotateRegions -bed regionToAnnotate.bed -ref referenceAnnotation.tsv -ref_sep "0 1 2" -ref_symbol "4 5"\n`, tail))
	}

	refPos := returnRefPos()
	symbolPos := returnSymbolPos()

	if FILENAMEOUT == "" {
		ext := path.Ext(BEDFILENAME.String())
		FILENAMEOUT = fmt.Sprintf("%s.annotated%s",
			BEDFILENAME[:len(BEDFILENAME)-len(ext)], ext)
	}

	if WRITEREF && WRITEINTERSECT {
		panic(fmt.Sprintf("Error! options -write_ref and -intersect cannot be TRUE together. Please chose one!\n"))
	}

	if !WRITEDIFF {
		utils.LoadRefCustomFileWithSymbol(REFBEDFILENAME, REFSEP, symbolPos, refPos)
	}

	utils.LoadPeaksCustom(REFBEDFILENAME, REFSEP, refPos)
	utils.CreatePeakIntervalTreeCustom(refPos, REFSEP)

	scanBedFileAndAddAnnotation(refPos)

	if REPLACEINPUT {
		utils.Check(os.Remove(string(BEDFILENAME)))
		utils.Check(os.Rename(FILENAMEOUT, string(BEDFILENAME)))
		fmt.Printf("File: %s edited\n", BEDFILENAME)
	} else {
		fmt.Printf("File: %s created\n", FILENAMEOUT)
	}
}

func returnRefPos() (refPos [3]int) {
	splitChar := " "

	if strings.Count(REFPOS, ",") > 0 {
		splitChar = ","
	}

	refPosSplit := strings.Split(REFPOS, splitChar)

	if len(refPosSplit) != 3 {
		panic(fmt.Sprintf("Error with ref_pos argument: %s should be an array of 3 ints (such as: 0 1 2) \n", REFPOS))
	}

	var err1, err2, err3 error

	refPos[0], err1 = strconv.Atoi(refPosSplit[0])
	refPos[1], err2 = strconv.Atoi(refPosSplit[1])
	refPos[2], err3 = strconv.Atoi(refPosSplit[2])

	if err1 != nil || err2 != nil || err3 != nil {
		panic(fmt.Sprintf("Error with ref_pos argument: %s should be an array of 3 ints (such as: 0 1 2) \n", REFPOS))
	}

	return refPos
}


func returnSymbolPos() (SymbolPos []int) {

	splitChar := " "

	if strings.Count(SYMBOLPOS, ",") > 0 {
		splitChar = ","
	}

	symbolPosSplit := strings.Split(SYMBOLPOS, splitChar)

	var err error

	SymbolPos = make([]int, len(symbolPosSplit))

	for pos, symbol := range symbolPosSplit {
		SymbolPos[pos], err = strconv.Atoi(symbol)

		if err != nil  {
			panic(fmt.Sprintf("Error with ref_pos argument: %s should be an array of ints \n", SYMBOLPOS))
	}
	}

	return SymbolPos
}

func scanBedFileAndAddAnnotation(refPos [3]int) {
	var intervals []interval.IntInterface
	var oneInterval interval.IntInterface
	var split, symbols []string
	var line, peakstr, symbol string
	var isInside, isUnique, isUniqueRef bool
	var start, end, count int
	var err error
	var inttree *interval.IntTree
	var buffer bytes.Buffer
	var intrange interval.IntRange
	var uniqueBed map[string]bool
	var intervalCenter int
	var centerDistance, minCenterDistance float64
	var oneIntervalID uintptr
	var peakIntervalTreeObject utils.PeakIntervalTreeObject
	var uniqueSymbolDict map[string]bool
	var writer io.WriteCloser

	scanner, file := BEDFILENAME.ReturnReader(0)
	defer utils.CloseFile(file)

	if STDOUT {
		writer = os.Stdout
	} else {
		writer = utils.ReturnWriter(FILENAMEOUT)
	}

	defer utils.CloseFile(writer)

	tStart := time.Now()

	if UNIQ {
		uniqueBed = make(map[string]bool)
	}


	if UNIQREF {
		peakIntervalTreeObject = utils.CreatePeakIntervalTreeObjectFromFile(
			BEDFILENAME)
	}

	for scanner.Scan() {
		line = scanner.Text()

		if line[0] == '#' {
			continue
		}

		split = strings.Split(line, "\t")

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		if inttree, isInside = utils.CHRINTERVALDICT[split[0]];!isInside {
			if !IGNOREUNANNOATED || WRITEDIFF{
				buffer.WriteString(line)
				buffer.WriteRune('\n')
				count++
			}
			continue
		}

		intervals = inttree.Get(
			utils.IntInterval{Start: start, End: end})

		switch {
		case len(intervals) == 0 :
			if !IGNOREUNANNOATED || WRITEDIFF{
				buffer.WriteString(line)
				buffer.WriteRune('\n')
				count++
			}
			if !isInside {
				continue
			}

		case  WRITEDIFF:
			continue
		}

		if UNIQSYMBOL {
			uniqueSymbolDict = make(map[string]bool)
		}

		intervalCenter, minCenterDistance = 0, -1
		isUnique = false
		isUniqueRef = true

		for _, oneInterval = range intervals {
			intrange = oneInterval.Range()
			oneIntervalID = oneInterval.ID()
			intervalCenter = intrange.Start - (intrange.End - intrange.Start) / 2
			centerDistance = math.Abs(float64((start - (end - start) / 2) - intervalCenter))

			if UNIQ {
				if minCenterDistance < 0 || centerDistance < minCenterDistance {
					minCenterDistance = centerDistance

					peakstr, symbols = returnPeakStrAndSymbol(
						line,
						oneIntervalID,
						intrange.Start,
						intrange.End,
						refPos)
				} else {
					continue
				}
			} else {
				peakstr, symbols = returnPeakStrAndSymbol(
					line,
					oneIntervalID,
					intrange.Start,
					intrange.End,
					refPos)
			}

			if UNIQSYMBOL {
				if uniqueSymbolDict[symbol] {
					continue
				}

				uniqueSymbolDict[symbol] = true
			}

			if UNIQREF {
				isUniqueRef = checkifUniqueRef(peakstr,
					oneIntervalID, refPos,
					&peakIntervalTreeObject)
			}

			if UNIQ {
				if uniqueBed[peakstr] {
					continue
				}

				uniqueBed[peakstr] = true
				isUnique = true

			} else  {

				if isUniqueRef {
					for _, symbol = range symbols {
						buffer.WriteString(peakstr)
						buffer.WriteRune('\t')
						buffer.WriteString(symbol)
						buffer.WriteRune('\n')
						count++
					}
				}
			}
		}

		if UNIQ && isUnique && isUniqueRef {
			for _, symbol = range symbols {
				buffer.WriteString(peakstr)
				buffer.WriteRune('\t')
				buffer.WriteString(symbol)
				buffer.WriteRune('\n')
				count++
			}
		}

		if count >= 5000 {
			_, err = writer.Write(buffer.Bytes())
			utils.Check(err)
			buffer.Reset()
		}

	}

	_, err = writer.Write(buffer.Bytes())
	utils.Check(err)
	buffer.Reset()

	tDiff := time.Since(tStart)
	fmt.Printf("done in time: %f s \n", tDiff.Seconds())
}


func returnPeakStrAndSymbol(line string, id uintptr, start, end int, refPos [3]int) (
	peakstr string, symbols []string) {

	var peak utils.Peak

	peakstr = utils.INTERVALMAPPING[id]

	peak.StringToPeakWithPos(peakstr, refPos)
	symbols = utils.PEAKSYMBOLDICT[peak]

	switch {
	case WRITEINTERSECT:
		peakstr = fmt.Sprintf("%s\t%d\t%d", peak.Slice[0],
			start,
			end)
	case WRITEREF:
		peakstr = peak.PeakToString()
	default:
		peakstr = line

	}

	return peakstr, symbols
}

func checkifUniqueRef(peakstr string, refpeakID uintptr, refPos [3]int,
	intervalObject *utils.PeakIntervalTreeObject ) bool {

	if !UNIQREF {
		return true
	}

	var refpeak utils.Peak
	var topID uintptr
	var centerDistance, minCenterDistance float64
	var intervalCenter, tss int
	var toppeakstr, refpeakstr string

	refpeakstr = utils.INTERVALMAPPING[refpeakID]

	refpeak.StringToPeakWithPos(refpeakstr, refPos)

	inttree := intervalObject.Chrintervaldict[refpeak.Chr()]

	intervals := inttree.Get(
		utils.IntInterval{Start: refpeak.Start, End: refpeak.End})

	intervalCenter, minCenterDistance = 0, -1

	tss = refpeak.Start + (refpeak.End - refpeak.Start) / 2

	for _, oneInterval := range intervals {
		intrange := oneInterval.Range()
		intervalCenter = intrange.Start - (intrange.End - intrange.Start) / 2

		centerDistance = math.Abs(float64((tss) - intervalCenter))

		if minCenterDistance < 0 || centerDistance < minCenterDistance {
			minCenterDistance = centerDistance
			topID = oneInterval.ID()
		}
	}

	toppeakstr = (*intervalObject).Intervalmapping[topID]

	return toppeakstr == peakstr
}
