package main

import (
	"log";
	"strings";
	"fmt";
	"strconv";
	"io"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"path"

)

/*INDEXTOOUTPUTNAME just the names of files */
var INDEXTOOUTPUTNAME map[string]map[string]string

func loadOutputFileIndex() {
	var line, index, indexstr, outtag string
	var split []string
	var isInside bool
	var outFilenameR1, outFilenameR2 string
	var r1r2Writers R1R2Writers

	if OUTPUTINDEXFILE == "" {
		return
	}

	INDEXTOOUTPUT = make(map[string]map[string]R1R2Writers)
	INDEXTOOUTPUTNAME = make(map[string]map[string]string)

	reader, file := utils.ReturnReader(OUTPUTINDEXFILE, 0)
	defer utils.CloseFile(file)

	ext := path.Ext(FASTQ_R1)


	for reader.Scan() {
		line = strings.Trim(reader.Text(), "\n\t\r")
		split = strings.Split(line, "\t")

		if len(split) != 3 {
			panic(fmt.Sprintf("Line: %s from file: %s should be splitted in 3",
				line, OUTPUTINDEXFILE))
		}

		index, indexstr, outtag = split[0], split[1], split[2]

		index = strings.Trim(index, " \t\n\r")
		indexstr = strings.Trim(indexstr, " \t\n\r")
		outtag = strings.Trim(outtag, " \t\n\r")

		if _, isInside = LENGTHDIC[index];!isInside {
			panic(fmt.Sprintf(
				"Error index: %s should be amongst %v", index, LENGTHDIC))
		}

		if _, isInside = INDEXTOOUTPUT[index];!isInside {
			INDEXTOOUTPUT[index] = make(map[string]R1R2Writers)
			INDEXTOOUTPUTNAME[index] = make(map[string]string)
		}

		if _, isInside = INDEXTOOUTPUT[index][indexstr];isInside {
			panic(
				fmt.Sprintf("Warning TAG: %s already exists for index: %s (targets: %s)\n",
				indexstr, index, INDEXTOOUTPUTNAME[index][indexstr]))
			continue
		}

		outFilenameR1 = fmt.Sprintf("%s/%s%s_R1.%s%s",
			OUTPUT_PATH, outtag, OUTPUT_TAG_NAME, OUTPUTFILETYPE, ext)
		outFilenameR2 = fmt.Sprintf("%s/%s%s_R2.%s%s",
			OUTPUT_PATH, outtag, OUTPUT_TAG_NAME, OUTPUTFILETYPE, ext)

		INDEXTOOUTPUTNAME[index][indexstr] = fmt.Sprintf(
			"%s/%s%s_R[12].%s%s",
			OUTPUT_PATH, outtag,
			OUTPUT_TAG_NAME, OUTPUTFILETYPE, ext)

		r1r2Writers = [2]io.WriteCloser{}

		r1r2Writers[0] = utils.ReturnWriter(outFilenameR1)
		r1r2Writers[1] = utils.ReturnWriter(outFilenameR2)

		INDEXTOOUTPUT[index][indexstr] = r1r2Writers
	}
}


func closeFileIndex() {
	for key := range INDEXTOOUTPUT {
		for barcode := range INDEXTOOUTPUT[key] {
			utils.CloseFile(INDEXTOOUTPUT[key][barcode][0])
			utils.CloseFile(INDEXTOOUTPUT[key][barcode][1])
		}
	}
}

func printSuccessFileIndex() {
	for key := range INDEXTOOUTPUT {
		for barcode := range INDEXTOOUTPUT[key] {
			fmt.Printf("Files %s created\n", INDEXTOOUTPUTNAME[key][barcode])
		}
	}
}

func loadIndexes(fnameList []string, dict * map[string]map[string]bool, reportName string) {
	countdict := make(map[string]int)

	(*dict) = make(map[string]map[string]bool)

	if len(fnameList) == 0 {
		fmt.Printf("#### index file not found! setting all tag length to: %d ####\n", TAGLENGTH)
		LENGTHDIC = map[string] int{
			"i5":TAGLENGTH,
			"i7":TAGLENGTH,
			"p5":TAGLENGTH,
			"p7":TAGLENGTH,
		}
		return
	}

	usedreadsfname := fmt.Sprintf("%s/used_barcodes%s.txt", OUTPUT_PATH, reportName)
	usedreadsf := utils.ReturnWriter(usedreadsfname)
	defer utils.CloseFile(usedreadsf)

	var tagNumber int
	var err error

	for _, fname := range(fnameList){
		scanner, file := utils.ReturnReader(fname, 0)
		defer utils.CloseFile(file)

	loop:
		for scanner.Scan() {
			line := scanner.Text()
			line = strings.Trim(line, "\n")
			split := strings.Split(line, "\t")

			if len(split) < 2 {
				panic(fmt.Sprintf("line %s in index: %s not conform!", line, fname))
			}

			if split[0] != "i5" && split[0] != "i7" &&  split[0] != "p5" && split[0] != "p7" {
				fmt.Printf("header not recognised in line %s skipping\n", split[0])
				continue loop
			}

			tagid, tagstring := split[0], split[1]
			length, isInside := LENGTHDIC[tagid]

			switch {
			case  !isInside:
				panic(fmt.Sprintf("tag ID %s is not valid (should be i5, p5, i7, p7)!", tagid))
			case len(tagstring) <= MAX_NB_MISTAKE_DICT[tagid]:
				panic(fmt.Sprintf("tag string %s not conform!", tagstring))
			case length == 0:
				LENGTHDIC[tagid] = len(tagstring)
			case length > 0 && len(tagstring) != length:
				panic(fmt.Sprintf("tag string %s for tag id %s has different" +
					" length than previous tags with similar id!", tagstring, tagid))
			}

			if _, isInside := (*dict)[tagid]; !isInside {
				(*dict)[tagid] = make(map[string]bool)
			}

			countdict[tagid]++

			if len(split) >= 3 {
				tagNumber, err = strconv.Atoi(split[2])
				if err != nil {
					tagNumber = countdict[tagid]
				}

			} else {
				tagNumber = countdict[tagid]
			}

			if len(INDEXESRANGE[tagid]) > 0 {
				if INDEXESRANGE[tagid][tagNumber] {
					(*dict)[tagid][tagstring] = true
					usedreadsf.Write(
						[]byte(fmt.Sprintf("%s\t%s\t%d\n", tagid, tagstring, tagNumber)))
				}

			} else {
				(*dict)[tagid][tagstring] = true
				usedreadsf.Write(
					[]byte(fmt.Sprintf("%s\t%s\t%d\n", tagid, tagstring, tagNumber)))
			}
		}
	}

}

/*LoadIndexRange load index range */
func LoadIndexRange() {
	INDEXESRANGE = make(map[string]map[int]bool)

	switch{
	case I5PLATES != "" && I5RANGE != "":
		log.Fatal("error I5_plates and I5_range both defined!")
	case I5PLATES != "":
		loadIndexRangeFrom("i5", I5PLATES, PLATESIZE)
	case I5RANGE != "":
		loadIndexRangeFrom("i5", I5RANGE, 1)
	}

	switch{
	case P7PLATES != "" && P7RANGE != "":
		log.Fatal("error P7_plates and P7_range both defined!")
	case P7PLATES != "":
		loadIndexRangeFrom("p7", P7PLATES, PLATESIZE)
	case P7RANGE != "":
		loadIndexRangeFrom("p7", P7RANGE, 1)
	}
}

func loadIndexRangeFrom(indexType string, str string, platesize int) {
	splitc := strings.Split(str, ",")
	INDEXESRANGE[indexType] = make(map[int]bool)

	for _ ,s := range splitc {
		if strings.Contains(s, "-") {
			ranges := strings.Split(s, "-")

			if len(ranges) != 2 {
				log.Fatal("Error i5_plates range wrong format!")
			}

			begin, err := strconv.Atoi(ranges[0])
			Check(err)
			end, err := strconv.Atoi(ranges[1])
			Check(err)

			for i:= begin; i <= end; i++ {
				for j := (i-1) * platesize; j < i * platesize; j++ {
					INDEXESRANGE[indexType][j+1] = true
				}
			}
		} else {

			index, err := strconv.Atoi(s)
			Check(err)
			for j := (index-1) * platesize; j < index * platesize; j++ {
				INDEXESRANGE[indexType][j+1] = true
			}
		}
	}
}

func checkIndexes(
	index_p7 string,
	index_i7 string,
	index_p5 string,
	index_i5 string) (
		bool,
		int,
		string,
		string,
		string,
		string) {

	if USENOINDEX{
		return true, 1, index_p7, index_i7, index_p5, index_i5
	}

	success_p7, ref_p7, replicate_p7 := checkOneIndex(index_p7, "p7")
	success_i7, ref_i7, replicate_i7 := checkOneIndex(index_i7, "i7")
	success_i5, ref_i5, replicate_i5 := checkOneIndex(index_i5, "i5")
	success_p5, ref_p5, replicate_p5 := checkOneIndex(index_p5, "p5")

	sum := replicate_p5 + replicate_i5 + replicate_i7 + replicate_p7
	success := (success_p7 && success_i7 && success_i5 && success_p5 && sum !=0)

	if !success{
		return false, 0, index_p7, index_i7, index_p5, index_i5
	}

	var repl int

	switch{
	case sum > 0:
		repl= 1
	default:
		repl= 2
	}

	return true, repl, ref_p7, ref_i7, ref_p5, ref_i5
}


func checkOneIndex(index string, indexName string) (success bool, ref string, replicate int) {

	switch{
	case index == "":
		return true, index, 0
	case USEALLBARCODES[indexName]:
		switch strings.Count(index, "N") {
		case 0:
			return true, index, 1
		default:
			goto fail
		}

	case INDEX_NO_REPLICATE:
		_, isInside := INDEX_NO_DICT[indexName][index]

		if isInside {
			return true, index, 1
		}

		isCloseEnougth, refindex := getNbMistake(index, indexName, INDEX_NO_DICT)

		if isCloseEnougth {
			return true, refindex, 1
		}

	default:
		_, isInside_r1 := INDEX_R1_DICT[indexName][index]
		_, isInside_r2 := INDEX_R2_DICT[indexName][index]

		switch{
		case isInside_r1 && isInside_r2:
			return true, index, 0
		case isInside_r1:
			return true, index, 1
		case isInside_r2:
			return true, index, -1
		default:
			isCloseEnougth_r1, refindex_r1 := getNbMistake(index, indexName, INDEX_R1_DICT)
			isCloseEnougth_r2, refindex_r2 := getNbMistake(index, indexName, INDEX_R2_DICT)

			switch{
			case isCloseEnougth_r1 && isCloseEnougth_r2:
				if refindex_r1 != refindex_r2 {
					goto fail
				}
				return true, refindex_r1, 0

			case isCloseEnougth_r1:
				return true, refindex_r1, 1
			case isCloseEnougth_r2:
				return true, refindex_r2, -1
			default:
				goto fail
			}
		}

	}
	fail:
	return false, index, 0
}

func getNbMistake (index string, indexName string, dict map[string]map[string]bool) (bool, string) {
	var nbMistake int

	for indexref := range dict[indexName] {
		nbMistake = 0
		for i, char := range index {
			if char != rune(indexref[i]) {
				nbMistake++
			}
		}

		if nbMistake <= MAX_NB_MISTAKE_DICT[indexName] {
			return true, indexref

		}
	}
	return false, ""
}
