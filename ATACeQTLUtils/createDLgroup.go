package main

import (
	"path"
	"fmt"
	"log"
	"time"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"strings"
	"strconv"
	"bytes"
	"sort"
)


func createDLGroup() {

	SNPDICT = make(map[string]map[int]float64)
	SNPIDTOSNPDICT = make(map[string]snpID)
	SNPTOSNPIDDICT = make(map[snpID]string)

	switch {

	case DBSNPFILE == "":
		log.Fatal(fmt.Printf("!!!! Error -dbSNP_file option must be provided\n"))
	case len(SNPFILES) == 0:
		log.Fatal(fmt.Printf("!!!! Error at least one -eQTL_gene_pair must be provided\n"))
	case OUTFILE == "":
		OUTFILE = fmt.Sprintf("%s.DL_group", DBSNPFILE)
	}

	ext := path.Ext(OUTFILE)

	OUTFILE = OUTFILE[:len(OUTFILE)-len(ext)]

	scanEQTLFile()
	scanDBSNPFile()
	writeUnknownSNP()
	scanDLFile()
}


func scanEQTLFile() {
	var line, chrID string
	var split, snv []string
	var pos, count int
	var err error
	var score float64
	var isInside bool

	tStart := time.Now()

	for _, fname := range SNPFILES{

		fmt.Printf("Scanning %s ...\n", fname)

		scanner, file := utils.ReturnReader(fname, 0)
		defer utils.CloseFile(file)
		scanner.Scan()

		for scanner.Scan() {
			line = scanner.Text()
			split = strings.Split(line, "\t")

			if len(split) < 12 {
				log.Fatal(fmt.Sprintf(
					"Error cannot split line: %s in more than 12 part using\t\n", split))
			}

			snv = strings.Split(split[0], "_")
			pos, err = strconv.Atoi(snv[1])
			utils.Check(err)

			chrID = snv[0]

			if _, isInside = SNPDICT[chrID];!isInside {
				SNPDICT[chrID] = make(map[int]float64)
			}

			score, err = strconv.ParseFloat(split[6], 64)
			utils.Check(err)

			SNPDICT[chrID][pos] = score
			count++
		}
	}

	fmt.Printf("Number of eQTL entries: %d\n", count)
	tDiff := time.Since(tStart)
	fmt.Printf("Scanning eQTL file done in time: %f s \n", tDiff.Seconds())
}


/*writeUnknownSNP Write SNP not found in dbSNP */
func writeUnknownSNP() {
	var isInside bool
	var buffer bytes.Buffer
	var snp snpID
	var count, pos int
	var err error

	fout := fmt.Sprintf("%s.unknownSNP.tsv", OUTFILE)
	tStart := time.Now()

	writer := utils.ReturnWriter(fout)

	for chrID := range SNPDICT {
		snp.chrID = chrID

		for pos = range SNPDICT[chrID] {
			snp.pos = pos

			if _, isInside = SNPTOSNPIDDICT[snp];!isInside {
				buffer.WriteString(snp.chrID)
				buffer.WriteRune('\t')
				buffer.WriteString(strconv.Itoa(snp.pos))
				buffer.WriteRune('\n')
				count++
			}
		}
		_, err = writer.Write(buffer.Bytes())
		utils.Check(err)
		buffer.Reset()
	}

	fmt.Printf("Number of unknown SNP %d\n", count)
	tDiff := time.Since(tStart)
	fmt.Printf("File: %s written in time: %f s \n", fout, tDiff.Seconds())
	fmt.Printf("Unknown SNP written in %s\n", fout)

}


func scanDLFile() {
	var line, snp1, snp2 string
	var split []string
	var isInside bool
	var snp snpID
	var score float64
	var snpItem item
	var count int

	snpScores := make([]item, 0, len(SNPIDTOSNPDICT))

	usedSnp := make(map[string]bool)
	snpGroup := make(map[string][]string)

	scanner, file := LDLINKFNAME.ReturnReader(0)
	defer utils.CloseFile(file)

	tStart := time.Now()

	for scanner.Scan() {
		line = scanner.Text()
		split = strings.Split(line, SEP)

		if len(split) < 2 {
			fmt.Printf("Error unable to split the line in two! line: %s\n", line)
			continue
		}

		snp1 = split[0]
		snp2 = split[1]

		if _, isInside = SNPIDTOSNPDICT[snp1]; !isInside {
			continue
		}

		if _, isInside = SNPIDTOSNPDICT[snp2]; !isInside {
			continue
		}

		snp = SNPIDTOSNPDICT[snp1]
		score = SNPDICT[snp.chrID][snp.pos]

		if _, isInside = usedSnp[snp1]; !isInside {
			usedSnp[snp1] = true

			snpItem.key = snp1
			snpItem.value = score

			snpScores = append(snpScores, snpItem)
		}

		count++
		snpGroup[snp1] = append(snpGroup[snp1], snp2)
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Scanning of LD file done in time: %f s \n", tDiff.Seconds())
	fmt.Printf("Number of links intersecting with used SNP %d\n", count)

	processScanDLResults(snpScores, &snpGroup)

}


func processScanDLResults(snpScores []item, snpGroup * map[string][]string) {

	var snp1, snp2 string
	var snpItem item
	var snp snpID
	var isInside bool
	var count1, count2 int
	var buffer1 bytes.Buffer
	var buffer2 bytes.Buffer
	var keys []string

	file1 := fmt.Sprintf("%s.dbSNPID.tsv", OUTFILE)
	file2 := fmt.Sprintf("%s.genomicID.tsv", OUTFILE)

	usedSnp := make(map[string]bool)

	fmt.Printf("Sorting...\n")
	tStart := time.Now()
	sort.Slice(snpScores, func(i, j int) bool {
		return snpScores[i].value > snpScores[j].value
	})
	tDiff := time.Since(tStart)
	fmt.Printf("Sorting done in time: %f s \n", tDiff.Seconds())

	for _, snpItem = range snpScores {
		snp1 = snpItem.key

		if _, isInside = usedSnp[snp1]; isInside {
			continue
		}

		for _, snp2 = range (*snpGroup)[snp1] {
			usedSnp[snp2] = true
		}

		usedSnp[snp1] = true
		keys = append(keys, snp1)
	}

	writer1 := utils.ReturnWriter(file1)
	defer utils.CloseFile(writer1)

	writer2 := utils.ReturnWriter(file2)
	defer utils.CloseFile(writer2)

	for _, snp1 = range keys {
		buffer1.WriteString(snp1)
		count1++
		count2++

		snp = SNPIDTOSNPDICT[snp1]
		buffer2.WriteString(snp.chrID)
		buffer2.WriteRune('_')
		buffer2.WriteString(strconv.Itoa(snp.pos))

		for _, snp2 = range (*snpGroup)[snp1] {
			buffer1.WriteRune('\t')
			buffer1.WriteString(snp2)

			snp = SNPIDTOSNPDICT[snp2]
			buffer2.WriteRune('\t')
			buffer2.WriteString(snp.chrID)
			buffer2.WriteRune('_')
			buffer2.WriteString(strconv.Itoa(snp.pos))

			count2++
		}

		buffer1.WriteRune('\n')
		buffer2.WriteRune('\n')

	}

	_, err := writer1.Write(buffer1.Bytes())
	utils.Check(err)

	_, err = writer2.Write(buffer2.Bytes())
	utils.Check(err)


	buffer1.Reset()
	buffer2.Reset()

	fmt.Printf("Number of DL group found: %d\n", count1)
	fmt.Printf("Number of SNP inside found: %d\n", count2)
	fmt.Printf("File %s written!\n", file1)
	fmt.Printf("File %s written!\n", file2)
}
