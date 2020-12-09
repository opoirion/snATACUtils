package main

import (
	"flag"
	utils "github.com/opoirion/snATACUtils/ATACdemultiplexUtils"
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

/*BEDPOS position of reference region in the -bed file */
var BEDPOS string

/*BEDPOSINT  bed pos int*/
var BEDPOSINT []int

/*SYMBOLPOS generic string or position of the coulmns used for annotations in the -ref file */
var SYMBOLPOS string

/*STDOUT write to stdout*/
var STDOUT bool

/*ANNOTATELINE annotate the full line rather than peak region*/
var ANNOTATELINE bool

/*SCOREFILTERCOLUMNS column to use to filter peaks*/
var SCOREFILTERCOLUMNS int

/*UNIQUEPEAKTOSYMBOL map used to link peak to unique top symbol */
var UNIQUEPEAKTOSYMBOL map[string]string

func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `
#################### MODULE TO ANNOTATE GENOMIC REGIONS FROM BED FILES ########################
This software presents some similarities with bedtools usage however it provides better customisations for bed file annotation when comparing two bed files with interesecting regions

"""Annotate bed file using a reference bed file containing the annotations"""
USAGE: ATACAnnotateRegions -bed <file> -ref <file> (optional -out <string> -unique -unique_ref -intersect -write_ref -edit -ref_sep "[3]int" -ref_symbol "[]int|str" -diff -stdout -annotate_line -score_pos)

for -ref_sep and -ref_symbol options, the input should be a string of numbers separated by whitespace and delimited with ". -ref_sep needs exactly three positions: 1) for the chromosomes column, 2) for the begining and 3) for the end of the region

Example: ATACAnnotateRegions -bed regionToAnnotate.bed -ref referenceAnnotation.tsv -ref_sep 0,1,2 -ref_symbol "4 5"

-ref_sep and -symbol_pos can be blank with " or comma separated: i.e."0 1 2" or 0,1,2. Also symbol_pos an be a generic string to annotate all ref regions

-score_pos is an alternative mechanism to retain unique peaks (no dupplicate) from the bed file by using an additional column of the bed file as score (needs to be float). if multiple entries exist for a given peak in the bed file, only the one with the highest score will be kept.

Here the three first columns of referenceAnnotation.tsv will be used to identify chromosome (column 0), start (column 1), and end (column 2) of each region, and regionToAnnotate.bed will be annotatd using columns 4 and 5 from referenceAnnotation.tsv
`)
		flag.PrintDefaults()
	}
	uniqsymbolstr := ""

	flag.Var(&BEDFILENAME, "bed", "name of the bed file no annotate")
	flag.Var(&REFBEDFILENAME, "ref", "name of the reference bed file containing the annotations")
	flag.StringVar(&FILENAMEOUT, "out", "", "name the output file(s)")
	flag.BoolVar(&REPLACEINPUT, "edit", false, `edit input bed file instead of creating a new file`)
	flag.BoolVar(&IGNOREUNANNOATED, "ignore", false, `ignore unnatotated peak`)
	flag.BoolVar(&WRITEINTERSECT, "intersect", false, `write intersection only`)
	flag.BoolVar(&UNIQ, "unique", false, `write only unique output peaks`)
	flag.BoolVar(&ANNOTATELINE, "annotate_line", false, `annotate the full line rather than the defined peak region`)
	flag.BoolVar(&UNIQREF, "unique_ref", false, `write only unique reference using the closest peak`)
	flag.IntVar(&SCOREFILTERCOLUMNS, "score_pos", -1, "(Require int) If used, refers to the column position containing score (float) to keep only the top unique link ")
	flag.StringVar(&uniqsymbolstr, "unique_symbols", "true", `write only unique symbols per peak`)
	flag.BoolVar(&WRITEDIFF, "diff", false, `write bed region if no intersection is found`)
	flag.BoolVar(&WRITEREF, "write_ref", false, `write bed region from reference file`)
	flag.BoolVar(&STDOUT, "stdout", false, `write to stdout`)
	flag.StringVar(&REFSEP, "ref_sep", "\t", "separator to define the bed region for the ref file")
	flag.StringVar(&REFPOS, "ref_pos", "", "separator to the bed region in ref for the ref file. Default: (0,1,2 for bed and 0,1,2,3,4,5 for bedpe files")
	flag.StringVar(&BEDPOS, "bed_pos", "", "separator to the bed region(s) in genomic coordinates for the bed file. Default: (0,1,2 for bed and 0,1,2,3,4,5 for bedpe files")
	flag.StringVar(&SYMBOLPOS, "symbol_pos", "3", "separator to the bed region in ref for the ref file")

	flag.Parse()

	if UNIQ && SCOREFILTERCOLUMNS > -1 {
		UNIQUEPEAKTOSYMBOL = make(map[string]string)

	} else {
		SCOREFILTERCOLUMNS = -1
	}

	UNIQSYMBOL = uniqsymbolstr == "true"

	tail := flag.Args()

	if len(tail) > 0 {
		panic(fmt.Sprintf(`Error wrongly formatted arguments: %s for -ref_sep and -ref_symbol options, the input should be a string of numbers separated by whitespace and delimited with ". -ref_sep needs exactly three positions: 1) for the chromosomes column, 2) for the begining and 3) for the end of the region

Example: ATACAnnotateRegions -bed regionToAnnotate.bed -ref referenceAnnotation.tsv -ref_sep "0 1 2" -ref_symbol "4 5"\n`, tail))
	}

	symbol := returnSymbolType()

	ext := path.Ext(BEDFILENAME.String())

	if BEDPOS == "" {
		switch ext {
		case ".bedpe":
			BEDPOS = "0,1,2,3,4,5"
		default:
			BEDPOS = "0,1,2"
		}
	}

	extRef := path.Ext(REFBEDFILENAME.String())

	if REFPOS == "" {
		switch extRef {
		case ".bedpe":
			REFPOS = "0,1,2,3,4,5"
		default:
			REFPOS = "0,1,2"
		}
	}

	refPos := returnPosIntSlice(REFPOS)

	if FILENAMEOUT == "" {
		FILENAMEOUT = fmt.Sprintf("%s.annotated%s",
			BEDFILENAME[:len(BEDFILENAME)-len(ext)], ext)
	}

	if WRITEREF && WRITEINTERSECT {
		panic(fmt.Sprintf("Error! options -write_ref and -intersect cannot be TRUE together. Please chose one!\n"))
	}

	if !WRITEDIFF {
		utils.LoadRefCustomFileWithSymbol(
			REFBEDFILENAME, REFSEP, symbol, refPos, SCOREFILTERCOLUMNS)
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

func returnPosIntSlice(refpos string) (refPos []int) {
	splitChar := " "

	if strings.Count(refpos, ",") > 0 {
		splitChar = ","
	}

	refPosSplit := strings.Split(refpos, splitChar)

	if len(refPosSplit) % 3 != 0 {
		panic(fmt.Sprintf(
			"Error with bed positional argument: refPos/bedPos: %s. Shouldb be an array with a length multiple of 3 ints separated by blank or , \n", refpos))
	}

	refPos = make([]int, len(refPosSplit))

	for i := 0; i < len(refPosSplit) / 3; i++ {
		j := i * 3
		var err1, err2, err3 error

		refPos[0 + j], err1 = strconv.Atoi(refPosSplit[0 + j])
		refPos[1 + j], err2 = strconv.Atoi(refPosSplit[1 + j])
		refPos[2 + j], err3 = strconv.Atoi(refPosSplit[2 + j])

		if err1 != nil || err2 != nil || err3 != nil {
			panic(fmt.Sprintf(
				"Error with bed positional argument: %s should be an array of ints (such as: 0 1 2) \n", refpos))
		}
	}

	return refPos
}


func writeDefault(line string, peak utils.Peak, buffer *bytes.Buffer) (count int) {

	if !ANNOTATELINE {
		line = peak.PeakToString()
	}

	buffer.WriteString(line)

	if WRITEDIFF {
		buffer.WriteRune('\n')
	} else {
		buffer.WriteString("\t\n")
	}
	count++

	return count
}

func returnSymbolType() (symbol utils.SymbolType) {

	splitChar := " "

	if strings.Count(SYMBOLPOS, ",") > 0 {
		splitChar = ","
	}

	symbolPosSplit := strings.Split(SYMBOLPOS, splitChar)

	var err error

	symbol.SymbolPos = make([]int, len(symbolPosSplit))

	for pos, sym := range symbolPosSplit {
		symbol.SymbolPos[pos], err = strconv.Atoi(sym)

		if err != nil  {
			symbol.SymbolPos = []int{}
			symbol.SymbolStr = SYMBOLPOS
			goto end
		}
	}

	end:
	return symbol
}

func scanBedFileAndAddAnnotation(refPosList []int) {
	var intervals []interval.IntInterface
	var oneInterval interval.IntInterface
	var symbols []string
	var refPos [3]int
	var line, peakstr string
	var isInside, isUnique, isUniqueRef bool
	var count int
	var err error
	var inttree *interval.IntTree
	var buffer bytes.Buffer
	var intrange interval.IntRange
	var uniqueBed map[string]bool
	var intervalCenter, nbPeak, nbPeakRef int
	var centerDistance, minCenterDistance float64
	var oneIntervalID uintptr
	var peakIntervalTreeObject utils.PeakIntervalTreeObject
	var writer io.WriteCloser
	var peak utils.Peak

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

	BEDPOSINT = returnPosIntSlice(BEDPOS)
	nbPeaksPerLine := utils.CheckIfPeakPosIsMutltipleOf3(BEDPOSINT)
	nbPeaksPerRef := utils.CheckIfPeakPosIsMutltipleOf3(refPosList)

	if UNIQREF {
		peakIntervalTreeObject = utils.CreatePeakIntervalTreeObjectFromFile(
			BEDFILENAME, "\t", BEDPOSINT)
	}

	for scanner.Scan() {
		line = scanner.Text()

		if line[0] == '#' {
			continue
		}

		for nbPeak = 0; nbPeak < nbPeaksPerLine; nbPeak++ {

			for nbPeakRef = 0; nbPeakRef < nbPeaksPerRef; nbPeakRef++ {
				refPos[0] = refPosList[0 + 3 * nbPeakRef]
				refPos[1] = refPosList[1 + 3 * nbPeakRef]
				refPos[2] = refPosList[2 + 3 * nbPeakRef]
				peak.StringToPeakWithPosAndStart(line, BEDPOSINT, nbPeak * 3)

				if inttree, isInside = utils.CHRINTERVALDICT[peak.Chr()];!isInside {
					if !IGNOREUNANNOATED || WRITEDIFF {
						count = writeDefault(line, peak, &buffer)
					}

					continue
				}

				intervals = inttree.Get(
					utils.IntInterval{Start: peak.Start, End: peak.End})

				switch {
				case len(intervals) == 0 :
					if !IGNOREUNANNOATED || WRITEDIFF {
						count = writeDefault(line, peak, &buffer)
					}

				case  WRITEDIFF:
					continue
				}

				intervalCenter, minCenterDistance = 0, -1
				isUnique = false
				isUniqueRef = true

				for _, oneInterval = range intervals {
					intrange = oneInterval.Range()
					oneIntervalID = oneInterval.ID()
					intervalCenter = intrange.Start - (intrange.End - intrange.Start) / 2
					centerDistance = math.Abs(float64((peak.Start - (peak.End - peak.Start) / 2) - intervalCenter))

					if UNIQ {
						if minCenterDistance < 0 || centerDistance < minCenterDistance {
							minCenterDistance = centerDistance

							peakstr, symbols = returnPeakStrAndSymbol(
								line,
								oneIntervalID,
								intrange.Start,
								intrange.End, nbPeak,
								refPos)
						} else {
							continue
						}
					} else {
						peakstr, symbols = returnPeakStrAndSymbol(
							line,
							oneIntervalID,
							intrange.Start,
							intrange.End, nbPeak,
							refPos)
					}

					if UNIQREF {
						isUniqueRef = checkifUniqueRef(
							peak.PeakToString(),
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
							if ANNOTATELINE {
								peakstr = line
							}

							count += writeToBuffer(
								symbols,
								peakstr,
								&buffer)
						}
					}
				}

				if !IGNOREUNANNOATED && !isUniqueRef {
					isUniqueRef = true
					symbols = []string{""}
				}

				if UNIQ && isUnique && isUniqueRef {
					if ANNOTATELINE {
						peakstr = line
					}

					count += writeToBuffer(
						symbols,
						peakstr,
						&buffer)
				}

				if count >= 5000 {
					_, err = writer.Write(buffer.Bytes())
					utils.Check(err)
					buffer.Reset()
				}

			}
		}

	}

	uniquePeakToSymbolToBuffer(&buffer, &writer)

	_, err = writer.Write(buffer.Bytes())
	utils.Check(err)
	buffer.Reset()

	tDiff := time.Since(tStart)

	if !STDOUT {
		fmt.Printf("done in time: %f s \n", tDiff.Seconds())
	}
}

func writeToBuffer(
	symbols []string,
	peakstr string,
	buffer * bytes.Buffer) (count int){
	var symbol string
	var uniqueSymbolDict map[string]bool

	if UNIQSYMBOL {
		uniqueSymbolDict = make(map[string]bool)
	}

	for _, symbol = range symbols {
		if UNIQSYMBOL {
			if uniqueSymbolDict[symbol] {
				continue
			}
			uniqueSymbolDict[symbol] = true
		}

		buffer.WriteString(peakstr)
		buffer.WriteRune('\t')
		buffer.WriteString(symbol)
		buffer.WriteRune('\n')
		count++

		if UNIQ {
			//only write first symbol per reference peaks
			break
		}
	}

	return count
}

func uniquePeakToSymbolToBuffer(buffer * bytes.Buffer, writer * io.WriteCloser) {
	count := 0
	var err error

	for peakstr, symbol := range UNIQUEPEAKTOSYMBOL {
		buffer.WriteString(peakstr)
		buffer.WriteRune('\t')
		buffer.WriteString(symbol)
		buffer.WriteRune('\n')
		count++

		if count > 5000 {
			_, err = (*writer).Write(buffer.Bytes())
			utils.Check(err)
			buffer.Reset()
			count = 0
		}
	}
}

func returnPeakStrAndSymbol(line string, id uintptr, start, end, nbPeak int, refPos [3]int) (
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
		peak.StringToPeakWithPosAndStart(line, BEDPOSINT, nbPeak * 3)
		peakstr = peak.PeakToString()
	}

	return peakstr, symbols
}

func checkifUniqueRef(peakstr string, refpeakID uintptr, refPos [3]int,
	intervalObject *utils.PeakIntervalTreeObject ) bool {

	if !UNIQREF {
		return true
	}

	var refpeak, toppeak utils.Peak
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
	toppeak.StringToPeak(toppeakstr)

	return peakstr == toppeak.PeakToString()
}
