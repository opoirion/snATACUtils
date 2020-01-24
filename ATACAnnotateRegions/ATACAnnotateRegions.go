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


func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, `
#################### MODULE TO ANNOTATE GENOMIC REGIONS FROM BED FILES ########################

"""Annotate bed file using a reference bed file containing the annotations"""
USAGE: ATACAnnotateRegions -chi2 -bed <fname> -peak <fname> -cluster <fname> (optionnal -out <string> -threads <int> -alpha <float> -write_all -split <int>)
`)
		flag.PrintDefaults()
	}

	flag.Var(&BEDFILENAME, "bed", "name of the bed file no annotate")
	flag.Var(&REFBEDFILENAME, "ref", "name of the reference bed file containing the annotations")
	flag.StringVar(&FILENAMEOUT, "out", "", "name the output file(s)")
	flag.BoolVar(&REPLACEINPUT, "edit", false, `edit input bed file instead of creating a new file`)
	flag.BoolVar(&IGNOREUNANNOATED, "ignore", false, `ignore unnatotated peak`)
	flag.BoolVar(&WRITEINTERSECT, "intersect", false, `write intersection only`)

	flag.Parse()

	if FILENAMEOUT == "" {
		ext := path.Ext(BEDFILENAME.String())
		FILENAMEOUT = fmt.Sprintf("%s.annotated%s",
			BEDFILENAME[:len(BEDFILENAME)-len(ext)], ext)
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
	var isInside bool
	var start, end, count int
	var err error
	var peak utils.Peak
	var inttree *interval.IntTree
	var buffer bytes.Buffer
	var intrange interval.IntRange

	scanner, file := BEDFILENAME.ReturnReader(0)
	defer utils.CloseFile(file)

	writer := utils.ReturnWriter(FILENAMEOUT)
	defer utils.CloseFile(writer)

	tStart := time.Now()

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

		for _, oneInterval = range intervals {

			peakstr = utils.INTERVALMAPPING[oneInterval.ID()]
			peak.StringToPeak(peakstr)

			symbol = utils.PEAKSYMBOLDICT[peak]

			if WRITEINTERSECT {
				intrange = oneInterval.Range()
				buffer.WriteString(peak[0])
				buffer.WriteRune('\t')
				buffer.WriteString(strconv.Itoa(intrange.Start))
				buffer.WriteRune('\t')
				buffer.WriteString(strconv.Itoa(intrange.End))
			} else {
				buffer.WriteString(line)
			}

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
