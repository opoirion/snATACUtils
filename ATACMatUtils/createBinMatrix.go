package main

import (
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"strings"
	"strconv"
	"fmt"
	"time"
	"bytes"
	"path"
)

type binPos struct {
	chr string
	index int
}

var BININDEX map[binPos]uint

/*BINSPARSEMATRIX map[cellID]map[bin]float64*/
var BINSPARSEMATRIX map[uint]map[uint]float64

/*BINSIZE bin size for bin matrix */
var BINSIZE int


func createBinSparseMatrix() {
	CELLIDCOUNT = make(map[string]int)

	loadCellIDDict(CELLSIDFNAME)
	initBinSparseMatrix()
	scanBedFileForBinMat()
	writeBinMatrixToFile(FILENAMEOUT)
}

func initBinSparseMatrix() {
	BINSPARSEMATRIX = make(map[uint]map[uint]float64)
	BININDEX = make(map[binPos]uint)

	for _, pos := range CELLIDDICT {
		BINSPARSEMATRIX[pos] = make(map[uint]float64)
	}
}

func scanBedFileForBinMat() {
	var line string
	var split []string
	var isInside bool
	var cellID, featureID, count uint
	var pos, index int
	var err error
	var bin binPos

	binList := []binPos{}
	tStart := time.Now()
	fmt.Printf("Scanning bed file...\n")

	bedReader, file := BEDFILENAME.ReturnReader(0)
	defer utils.CloseFile(file)

	for bedReader.Scan() {
		line = bedReader.Text()
		split = strings.Split(line, "\t")

		if cellID, isInside = CELLIDDICT[split[3]];!isInside {
			continue
		}

		pos, err = strconv.Atoi(split[1])
		utils.Check(err)

		index = (pos) / 5000

		bin.chr = split[0]
		bin.index = index

		if featureID, isInside = BININDEX[bin];!isInside {
			featureID = count
			BININDEX[bin] = count
			binList = append(binList, bin)
			count++
		}

		BINSPARSEMATRIX[cellID][featureID]++
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Scanning done in time: %f s \n", tDiff.Seconds())

	writeBinList(binList)
}


func writeBinList(binList []binPos) {
	var bin binPos
	var index int
	ext := path.Ext(CELLSIDFNAME.String())

	var buffer bytes.Buffer

	outfname := fmt.Sprintf("%s.ygi",
		CELLSIDFNAME[:len(CELLSIDFNAME) - len(ext)])

	writer := utils.ReturnWriter(outfname)
	defer utils.CloseFile(writer)

	for _, bin = range binList {
		index = int(bin.index) * bin.index

		buffer.WriteString(bin.chr)
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.Itoa(index))
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.Itoa(index + BINSIZE))
		buffer.WriteRune('\n')
	}

	fmt.Printf("File %s written!\n", outfname)
}


func writeBinMatrixToFile(outfile string) {
	fmt.Printf("writing Bin matrix to output file...\n")
	tStart := time.Now()

	var buffer bytes.Buffer
	var cellPos uint
	var index uint
	var binValue float64

	writer := utils.ReturnWriter(outfile)

	defer utils.CloseFile(writer)

	for cellPos = range BINSPARSEMATRIX {
		for index = range BINSPARSEMATRIX[cellPos] {
			buffer.WriteString(strconv.Itoa(int(cellPos)))
			buffer.WriteString(SEP)
			buffer.WriteString(strconv.Itoa(int(index)))
			buffer.WriteString(SEP)

			binValue = BINSPARSEMATRIX[cellPos][index]
			buffer.WriteString(fmt.Sprintf("%f\n", binValue))

			writer.Write(buffer.Bytes())
			buffer.Reset()
		}
	}

	fmt.Printf("file: %s created!\n", outfile)
	tDiff := time.Since(tStart)
	fmt.Printf("done in time: %f s \n", tDiff.Seconds())
}
