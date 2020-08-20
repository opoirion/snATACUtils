package main

import (
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"strings"
	"strconv"
	"fmt"
	"time"
	"bytes"
	"path"
	"sync"
	"github.com/biogo/store/interval"
	"sort"
)

type binPos struct {
	chr string
	index int
}

/*BININDEX dict for bin index: map[bin]index */
var BININDEX map[binPos]uint

/*BINSIZE bin size for bin matrix */
var BINSIZE int

/*BININDEXCOUNT index used to create  */
var BININDEXCOUNT uint

/*BININDEXMUTEX bin index mutex */
var BININDEXMUTEX sync.Mutex


func createBinSparseMatrix() {
	loadCellIDDict(CELLSIDFNAME)

	XGIDIM = len(CELLIDDICT)
	initIntSparseMatrix()
	BININDEX = make(map[binPos]uint)

	if PEAKFILE != "" {
		YGIDIM = utils.LoadPeaks(PEAKFILE)
		utils.CreatePeakIntervalTree()
		utils.InitIntervalDictsThreading(THREADNB)
		createBinSparseMatrixOneFileThreading(BEDFILENAME)
	} else {
		scanBedFileForBinMat()
	}


	switch MATRIXFORMAT {
	case mtx:
		writeIntMatrixToCOOFile(FILENAMEOUT, true)
	case coo:
		writeIntMatrixToCOOFile(FILENAMEOUT, false)
	case taiji:
		writeIntMatrixToTaijiFile(FILENAMEOUT, true)
	case dense:
		writeIntMatrixToDenseFile(FILENAMEOUT, true)
	case denseTranspose:
		writeIntMatrixToDenseTransposeFile(FILENAMEOUT)
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

		index = (pos) / BINSIZE

		bin.chr = split[0]
		bin.index = index

		if featureID, isInside = BININDEX[bin];!isInside {
			featureID = count
			BININDEX[bin] = count
			binList = append(binList, bin)
			count++
		}

		if NORM {
			TOTALREADSCELL[cellID]++
		}

		INTSPARSEMATRIX[cellID][featureID]++
	}

	YGIDIM = len(binList)

	tDiff := time.Since(tStart)
	fmt.Printf("Scanning done in time: %f s \n", tDiff.Seconds())

	writeBinList(binList)
}


func writeBinList(binList []binPos) {
	var bin binPos
	var index int
	var outfname string

	if YGIOUT != "" {
		outfname = YGIOUT
	} else {
		ext := path.Ext(CELLSIDFNAME.String())
		outfname = fmt.Sprintf("%s.ygi",
		CELLSIDFNAME[:len(CELLSIDFNAME) - len(ext)])
	}

	var buffer bytes.Buffer

	writer := utils.ReturnWriter(outfname)
	defer utils.CloseFile(writer)

	for _, bin = range binList {
		index = int(bin.index) * BINSIZE

		buffer.WriteString(bin.chr)
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.Itoa(index))
		buffer.WriteRune('\t')
		buffer.WriteString(strconv.Itoa(index + BINSIZE))
		buffer.WriteRune('\n')
	}

	_, err := writer.Write(buffer.Bytes())
	utils.Check(err)

	fmt.Printf("File %s written!\n", outfname)
}

func createBinSparseMatrixOneFileThreading(bedfilename utils.Filename) {
	var nbReads uint
	var bufferLine1 [BUFFERSIZE]string
	var bufferLine2 [BUFFERSIZE]string
	var bufferPointer * [BUFFERSIZE]string

	isBuffer1 := true
	bufferPointer = &bufferLine1

	var bufferIt int
	var waiting sync.WaitGroup

	bedReader, file := bedfilename.ReturnReader(0)

	defer utils.CloseFile(file)

	scanBed:
	for bedReader.Scan() {
		bufferPointer[bufferIt] = bedReader.Text()
		nbReads++
		bufferIt++


		if bufferIt >= BUFFERSIZE {
			chunk := bufferIt / THREADNB
			bufferStart := 0
			bufferStop := chunk

			for i := 0; i < THREADNB;i++{
				go updateBinSparseMatrixOneThread(bufferPointer , bufferStart, bufferStop, i, &waiting)

				bufferStart += chunk
				bufferStop += chunk

				if i == THREADNB - 1 {
					bufferStop = bufferIt
				}
			}

			bufferIt = 0

			if isBuffer1 {
				bufferPointer = &bufferLine2
				isBuffer1 = false
				goto scanBed
			} else {
				bufferPointer = &bufferLine1
				isBuffer1 = true
				waiting.Wait()
			}
		}
	}

	waiting.Wait()

	if bufferIt > 0 {
		waiting.Add(1)
		updateBinSparseMatrixOneThread(bufferPointer , 0, bufferIt, 0, &waiting)
	}

	sortBinIndexAndwrite()

}

func sortBinIndexAndwrite() {
	type binSort struct {
		bin binPos
		indexMat uint
	}

	binList := []binSort{}

	for bin, indexMat := range BININDEX {
		binList = append(binList, binSort{bin:bin, indexMat:indexMat})
	}

	sort.Slice(binList, func(i, j int) bool {
		return binList[i].indexMat < binList[j].indexMat
	})

	binListReady := make([]binPos, len(binList))

	for i, binsort := range binList {
		binListReady[i] = binsort.bin
	}

	writeBinList(binListReady)

}

func updateBinSparseMatrixOneThread(bufferLine * [BUFFERSIZE]string, bufferStart ,bufferStop, threadnb int,
	waiting * sync.WaitGroup) {
	defer waiting.Done()
	var split []string
	var isInside bool
	var start, end int
	var err error
	var bin binPos
	var index int
	var featureID, cellID uint

	waiting.Add(1)
	var intervals []interval.IntInterface

	type tmpMatUnit struct {
		bin binPos
		cellID uint
	}

	var tmpResult []tmpMatUnit
	var oneFeat tmpMatUnit

	for i := bufferStart; i < bufferStop;i++ {

		split = strings.Split(bufferLine[i], "\t")

		if cellID, isInside = CELLIDDICT[split[3]];!isInside {
			continue
		}

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		if _, isInside = utils.CHRINTERVALDICT[split[0]];!isInside {
			continue
		}

		intervals = utils.CHRINTERVALDICTTHREAD[threadnb][split[0]].Get(
			utils.IntInterval{Start: start, End: end})

		if len(intervals) == 0 {
			continue
		}

		index = (start) / BINSIZE

		bin.chr = split[0]
		bin.index = index

		oneFeat.bin = bin
		oneFeat.cellID = cellID

		tmpResult = append(tmpResult, oneFeat)
	}

	BININDEXMUTEX.Lock()

	for _, oneFeat = range tmpResult {
		if featureID, isInside = BININDEX[bin];!isInside {
			featureID = BININDEXCOUNT
			BININDEX[bin] = BININDEXCOUNT
			BININDEXCOUNT++
		}

		if NORM {
			TOTALREADSCELL[cellID]++
		}

		INTSPARSEMATRIX[cellID][featureID]++
	}

	BININDEXMUTEX.Unlock()
}
