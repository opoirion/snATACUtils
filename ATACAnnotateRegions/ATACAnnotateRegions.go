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

/*WRITEREF write_ref write bed region from reference file */
var WRITEREF bool


func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `
#################### MODULE TO ANNOTATE GENOMIC REGIONS FROM BED FILES ########################

"""Annotate bed file using a reference bed file containing the annotations"""
USAGE: ATACAnnotateRegions -bed <file> -ref <file> (optionnal -out <string> -unique -unique_ref -intersect -write_ref)
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
	flag.BoolVar(&WRITEREF, "write_ref", false, `write bed region from reference file`)

	flag.Parse()

	if FILENAMEOUT == "" {
		ext := path.Ext(BEDFILENAME.String())
		FILENAMEOUT = fmt.Sprintf("%s.annotated%s",
			BEDFILENAME[:len(BEDFILENAME)-len(ext)], ext)
	}

	if WRITEREF && WRITEINTERSECT {
		panic(fmt.Sprintf("Error! options -write_ref and -intersect cannot be TRUE together. Please chose one!\n"))
	}

	utils.LoadRefBedFileWithSymbol(REFBEDFILENAME)
	utils.LoadPeaks(REFBEDFILENAME)
	utils.CreatePeakIntervalTree()

	scanBedFileAndAddAnnotation()

	if REPLACEINPUT {
		utils.Check(os.Remove(string(BEDFILENAME)))
		utils.Check(os.Rename(FILENAMEOUT, string(BEDFILENAME)))
		fmt.Printf("File: %s edited\n", BEDFILENAME)
	} else {
		fmt.Printf("File: %s created\n", FILENAMEOUT)
	}
}


func scanBedFileAndAddAnnotation() {
	var intervals []interval.IntInterface
	var oneInterval interval.IntInterface
	var split []string
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

	scanner, file := BEDFILENAME.ReturnReader(0)
	defer utils.CloseFile(file)

	writer := utils.ReturnWriter(FILENAMEOUT)
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
		split = strings.Split(line, "\t")

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		if inttree, isInside = utils.CHRINTERVALDICT[split[0]];!isInside {
			continue
		}

		intervals = inttree.Get(
			utils.IntInterval{Start: start, End: end})

		if len(intervals) == 0 && !IGNOREUNANNOATED {
			buffer.WriteString(line)
			buffer.WriteRune('\n')
			count++
		}

		intervalCenter, minCenterDistance = 0, -1
		isUnique = false
		isUniqueRef = true

		for _, oneInterval = range intervals {
			intrange = oneInterval.Range()
			oneIntervalID = oneInterval.ID()
			intervalCenter = intrange.End - intrange.Start
			centerDistance = math.Abs(float64((start - end) - intervalCenter))

			if UNIQ {
				if minCenterDistance < 0 || centerDistance < minCenterDistance {
					minCenterDistance = centerDistance

					peakstr, symbol = returnPeakStrAndSymbol(
						line,
						oneIntervalID,
						intrange.Start,
						intrange.End)
				} else {
					continue
				}
			} else {
				peakstr, symbol = returnPeakStrAndSymbol(
					line,
					oneIntervalID,
					intrange.Start,
					intrange.End)
			}

			if UNIQREF {
				isUniqueRef = checkifUniqueRef(peakstr,
					oneIntervalID,
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
					buffer.WriteString(peakstr)
					buffer.WriteRune('\t')
					buffer.WriteString(symbol)
					buffer.WriteRune('\n')
					count++
				}
			}
		}

		if UNIQ && isUnique && isUniqueRef {
			buffer.WriteString(peakstr)
			buffer.WriteRune('\t')
			buffer.WriteString(symbol)
			buffer.WriteRune('\n')
			count++
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


func returnPeakStrAndSymbol(line string, id uintptr, start, end int) (
	peakstr, symbol string) {

	var peak utils.Peak

	peakstr = utils.INTERVALMAPPING[id]
	peak.StringToPeak(peakstr)
	symbol = utils.PEAKSYMBOLDICT[peak]

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

	return peakstr, symbol
}

func checkifUniqueRef(peakstr string, refpeakID uintptr,
	intervalObject *utils.PeakIntervalTreeObject ) bool {

	if !UNIQREF {
		return true
	}

	var refpeak utils.Peak
	var topID uintptr
	var centerDistance, minCenterDistance float64
	var intervalCenter int
	var toppeakstr, refpeakstr string

	refpeakstr = utils.INTERVALMAPPING[refpeakID]

	refpeak.StringToPeak(refpeakstr)

	inttree := intervalObject.Chrintervaldict[refpeak.Chr()]

	intervals := inttree.Get(
		utils.IntInterval{Start: refpeak.Start, End: refpeak.End})

	intervalCenter, minCenterDistance = 0, -1

	for _, oneInterval := range intervals {
		intrange := oneInterval.Range()
		intervalCenter = intrange.End - intrange.Start
		centerDistance = math.Abs(float64((refpeak.Start - refpeak.End) - intervalCenter))

		if minCenterDistance < 0 || centerDistance < minCenterDistance {
			minCenterDistance = centerDistance
			topID = oneInterval.ID()
		}
	}

	toppeakstr = (*intervalObject).Intervalmapping[topID]

	if toppeakstr == peakstr {
		return true
	}

	return false
}
