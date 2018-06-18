package main

import (
	"bufio"
	"os";
	"log";
	"strings";
	"fmt"
)


func loadIndexes(fname string, dict * map[string]map[string]bool) {
	if fname == "" {
		return
	}

	(*dict) = make(map[string]map[string]bool)

	file_open, err := os.OpenFile(fname, 0, 0)

	if err != nil {
		log.Fatal(err)
	}

	reader := bufio.NewReader(file_open)
	scanner := bufio.NewScanner(reader)

	first := true


	loop:
	for scanner.Scan() {
		line := scanner.Text()
		line = strings.Trim(line, "\n")
		split := strings.Split(line, "\t")

		if len(split) != 2 {
			panic(fmt.Sprintf("line %s in index: %s not conform!", line, fname))
		}

		if split[0] != "i5" && split[0] != "i7" &&  split[0] != "p5" && split[0] != "p7" {
			fmt.Printf("header not recognised in line %s skipping\n", split[0])
			continue loop
		}

		switch {
		case first && TAGLENGTH != len(split[1]):
			fmt.Printf(`taglength: %d different from the length of first library loaded: %s.
 Automaticcaly change the length\n`, TAGLENGTH, split[1])
		case !first && TAGLENGTH != len(split[1]):
			log.Fatal(fmt.Sprintf("error! different taglength found!: %s", split[1]))
		case first:
			first = false
		}

		tagid, tagstring := split[0], split[1]

		if _, isInside := (*dict)[tagid]; !isInside {
			(*dict)[tagid] = make(map[string]bool)
		}

		(*dict)[tagid][tagstring] = true

	}

}

func checkIndexes(
	index_p7 * string,
	index_i7 * string,
	index_p5 * string,
	index_i5 * string) (bool, int) {
	if INDEX_REPLICATE_R1 == "" && INDEX_NO_REPLICATE == ""{
		return true, 1
	}

	success_p7, toChange_p7, ref_p7, replicate_p7 := checkOneIndex(*index_p7, "p7")
	success_i7, toChange_i7, ref_i7, replicate_i7 := checkOneIndex(*index_i7, "i7")
	success_i5, toChange_i5, ref_i5, replicate_i5 := checkOneIndex(*index_i5, "i5")
	success_p5, toChange_p5, ref_p5, replicate_p5 := checkOneIndex(*index_p5, "p5")

	sum := replicate_p5 + replicate_i5 + replicate_i7 + replicate_p7
	success := (success_p7 && success_i7 && success_i5 && success_p5 && sum !=0)

	if !success{
		return false, 0
	}

	if toChange_p7{
		*index_p7 = ref_p7
	}

	if toChange_i7{
		*index_i7 = ref_i7
	}

	if toChange_p5{
		*index_p5 = ref_p5
	}

	if toChange_i5{
		*index_i5 = ref_i5
	}

	var repl int

	switch{
	case sum > 0:
		repl= 1
	default:
		repl= 2
	}

	return true, repl
}


func checkOneIndex(index string, indexName string) (bool, bool, string, int) {

	switch{
	case INDEX_NO_REPLICATE != "":
		_, isInside := INDEX_NO_DICT[indexName][index]

		if isInside {
			return true, false, "", 1
		}

		isCloseEnougth, refindex := getNbMistake(index, indexName, INDEX_NO_DICT)

		if isCloseEnougth {
			return true, true, refindex, 1
		}

	default:
		_, isInside_r1 := INDEX_R1_DICT[indexName][index]
		_, isInside_r2 := INDEX_R2_DICT[indexName][index]

		switch{
		case isInside_r1 && isInside_r2:
			return true, false, "", 0
		case isInside_r1:
			return true, false, "", 1
		case isInside_r2:
			return true, false, "", -1
		default:
			isCloseEnougth_r1, refindex_r1 := getNbMistake(index, indexName, INDEX_R1_DICT)
			isCloseEnougth_r2, refindex_r2 := getNbMistake(index, indexName, INDEX_R2_DICT)

			switch{
			case isCloseEnougth_r1 && isCloseEnougth_r2:
				if refindex_r1 != refindex_r2 {
					goto fail
				}
				return true, true, refindex_r1, 0

			case isCloseEnougth_r1:
				return true, true, refindex_r1, 1
			case isCloseEnougth_r2:
				return true, true, refindex_r2, -1
			default:
				goto fail
			}
		}

	}
	fail:
	return false, false, "", 0
}

func getNbMistake (index string, indexName string, dict map[string]map[string]bool) (bool, string) {
	var nbMistake int

	for indexref, _ := range dict[indexName] {
		nbMistake = 0
		for i, char := range index {
			if char != rune(indexref[i]) {
				nbMistake++
			}
		}

		if nbMistake <= MAX_NB_MISTAKE {
			return true, indexref

		}
	}
	return false, ""
}
