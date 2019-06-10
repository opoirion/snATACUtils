package main


import (
	"flag"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"github.com/biogo/store/interval"
	"fmt"
	"log"
	"strings"
	"strconv"
	"bytes"
	"path"
	"io"
	"sort"
	"time"
)


/*BEDFILENAME bed file name (input) */
var BEDFILENAME string

/*PEAKFILE file containg cluster-specific peaks */
var PEAKFILE string

/*OUTFILE output file */
var OUTFILE string

/*GENEFILE Specific gene file */
var GENEFILE string

/*SNPFILES multiple input files for GWAS SNP-genes */
var SNPFILES utils.ArrayFlags

/*GENEIDTONAMEFILE gene ID to name conversion File */
var GENEIDTONAMEFILE string

/*WRITESNPTOBED write significant SNP to bed file */
var WRITESNPTOBED bool

/*CHRINTERVALDICT chr ID <-> interval tree */
var CHRINTERVALDICT map[string]*interval.IntTree

/*PEAKCHRINTERVALDICT chr ID <-> interval tree */
var PEAKCHRINTERVALDICT map[string]*interval.IntTree

/*GENEDICT genename <string> -> score <float> */
var GENEDICT map[string]float64

/*GENEUINTDICT genename <string> -> score <float> */
var GENEUINTDICT map[uintptr]string

/*PEAKUNINTDICT genename <string> -> score <float> */
var PEAKUNINTDICT map[uintptr]string

/*CELLIDEQTLDICT genename <string> -> score <float> */
var CELLIDEQTLDICT map[string]map[string]bool

/*CELLIDDICT cell ID<->pos */
var CELLIDDICT map[string]bool

/*CELLIDCOMPLDICT cell ID<->pos */
var CELLIDCOMPLDICT map[string]bool

/*EGENEDICT cell ID<->pos */
var EGENEDICT map[string]string

/*NBCELLS number of cells in the clusters */
var NBCELLS int

/*NBCELLSCOMPL cell */
var NBCELLSCOMPL int

/*CELLIDFNAME cell ID<->pos */
var CELLIDFNAME string

/*CELLIDCOMPLFNAME cell ID<->pos */
var CELLIDCOMPLFNAME string

/*SNPPADDING Padding for SNP inclusion */
var SNPPADDING int

/*IntInterval Integer-specific intervals used for intervalTree*/
type IntInterval struct {
	Start, End int
	UID        uintptr
	Payload    interface{}
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


func main() {
	flag.BoolVar(&WRITESNPTOBED, "write_snp", false, "Write potential SNPs to a bed file")
	flag.IntVar(&SNPPADDING, "padding", 0, "Padding added before and after SNP location for SNP inclusion")
	flag.StringVar(&BEDFILENAME, "bed", "", "name of the bed file")
	flag.StringVar(&PEAKFILE, "peak_file", "", "Cluster-specific Peak file ")
	flag.StringVar(&CELLIDFNAME, "xgi", "", "cluster specific cell barcodes ")
	flag.StringVar(&CELLIDCOMPLFNAME, "xgi_compl", "", "All cell barcodes ")
	flag.StringVar(&GENEFILE, "gene_file", "", "Cluster-specific Gene file ")
	flag.StringVar(&GENEIDTONAMEFILE, "gene_ID_to_name", "", "File used for gene name conversion (gene_id\tgene_name) ")
	flag.Var(&SNPFILES, "GWAS_gene_pair", "name of one or multiple GWAS files (-GWAS_gene_pair <fname1> -GWAS_gene_pair <fname2>)")
	flag.StringVar(&OUTFILE, "out", "", "Name of the output file")

	flag.Parse()

	switch {

	case BEDFILENAME == "":
		log.Fatal(fmt.Printf("!!!! Error -bed option must be provided\n"))
	case PEAKFILE == "":
		log.Fatal(fmt.Printf("!!!! Error -peak option must be provided\n"))
	case GENEFILE == "":
		log.Fatal(fmt.Printf("!!!! Error -gene_file option must be provided\n"))
	case CELLIDFNAME == "":
		log.Fatal(fmt.Printf("!!!! Error -xgi must be provided\n"))
	case CELLIDCOMPLFNAME == "":
		log.Fatal(fmt.Printf("!!!! Error -xgi_compl must be provided\n"))
	case len(SNPFILES) == 0:
		log.Fatal(fmt.Printf("!!!! Error at least one -GWAS_gene_pair must be provided\n"))
	case OUTFILE == "":
		OUTFILE = fmt.Sprintf("%s.eQTL.tsv", BEDFILENAME)

	}

	tStart := time.Now()
	//Create and init cell ID dicts
	loadCellIDDicts()
	//Create Ref Gene Array
	createGeneArray()
	//Create gene ID <-> gene name dict
	createRefGeneDict()
	//Create peak interval tree
	createPeakIntervalTree()
	//Create SNP interval tree
	createSNPIntervalTree()
	//Scan the bed files
	scanBed()
	tDiff := time.Since(tStart)
	fmt.Printf("done in time: %f s \n", tDiff.Seconds())

}


func loadCellIDDicts() {
	CELLIDEQTLDICT = make(map[string]map[string]bool)

	CELLIDDICT = utils.LoadCellIDDict(CELLIDFNAME)
	NBCELLS = len(CELLIDDICT)

	CELLIDCOMPLDICT = utils.LoadCellIDDict(CELLIDCOMPLFNAME)

	for key := range CELLIDCOMPLDICT {
		if _, isInside := CELLIDEQTLDICT[key];!isInside {
			CELLIDEQTLDICT[key] = make(map[string]bool)
		}

		if CELLIDDICT[key] {
			delete(CELLIDCOMPLDICT, key)
		}
	}

	NBCELLSCOMPL = len(CELLIDCOMPLDICT)

	for key := range CELLIDDICT {
		if _, isInside := CELLIDEQTLDICT[key];!isInside {
			CELLIDEQTLDICT[key] = make(map[string]bool)
		}
	}
}

func scanBed() {
	var err error
	var buffer bytes.Buffer
	var keys []string
	var scoreNCl, scoreNCompl float64

	countDicts := scanOneBedFile(BEDFILENAME)
	writer := utils.ReturnWriter(OUTFILE)
	defer utils.CloseFile(writer)
	defer fmt.Printf("file: %s written\n", OUTFILE)

	for key := range countDicts {
		keys = append(keys, key)
	}

	sort.Strings(keys)

	writer.Write([]byte("#"))
	writer.Write([]byte("eQTL_linkage\t"))
	_, err = writer.Write([]byte(strings.Join(keys, "\t")))
	utils.Check(err)
	writer.Write([]byte("Score Normed Cluster\t"))
	writer.Write([]byte("Score Normed Compl.\t"))
	writer.Write([]byte("Nb. cell Cluster\t"))
	writer.Write([]byte("Nb. cell Compl.\t\n"))

	firstKey := keys[0]
	count := 0

	for key := range countDicts[firstKey] {
		buffer.WriteString(key)

		for _, mainKey := range keys {
			buffer.WriteRune('\t')
			buffer.WriteString(strconv.Itoa(countDicts[mainKey][key]))
		}

 		scoreNCl = float64(countDicts["Score cluster"][key]) / float64(NBCELLS)
		scoreNCompl = float64(countDicts["Score compl."][key]) / float64(NBCELLSCOMPL)

		buffer.WriteRune('\t')
		buffer.WriteString(fmt.Sprintf("%f", scoreNCl))
		buffer.WriteRune('\t')
		buffer.WriteString(fmt.Sprintf("%f", scoreNCompl))
		buffer.WriteRune('\t')
		buffer.WriteString(fmt.Sprintf("%d", NBCELLS))
		buffer.WriteRune('\t')
		buffer.WriteString(fmt.Sprintf("%d", NBCELLSCOMPL))
		buffer.WriteRune('\n')

		_, err = writer.Write(buffer.Bytes())
		utils.Check(err)
		buffer.Reset()
		count++
	}

	fmt.Printf("Number of eQTL enhancer-genes detected: %d\n", count)
}


func scanOneBedFile(bedfile string) map[string]map[string]int  {
	var line, chrID, cellID string
	var start, end int
	var inter IntInterval
	var split []string
	var err error
	var useGeneCountCompl, isInside, isInsideCompl bool

	countDicts := make(map[string]map[string]int)

	countDicts["Score cluster"] = make(map[string]int)
	countDicts["Score compl."] = make(map[string]int)
	countDicts["Cell number cluster"] = make(map[string]int)
	countDicts["Cell number compl."] = make(map[string]int)

	bedReader, file := utils.ReturnReader(BEDFILENAME, 0)
	defer utils.CloseFile(file)

	bedLoop:
	for bedReader.Scan() {
		line = bedReader.Text()
		split = strings.Split(line, "\t")

		cellID = split[3]

		_, isInside = CELLIDDICT[cellID]
		_, isInsideCompl = CELLIDCOMPLDICT[cellID]

		switch {
		case !isInside && !isInsideCompl:
			continue bedLoop
		case isInside:
			useGeneCountCompl = false
		case isInsideCompl:
			useGeneCountCompl = true
		case isInside && isInsideCompl:
			panic("Error cells inside cluster and compl!")
		}

		chrID = split[0][3:]

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)
		inter.Start = start
		inter.End = end

		if _, isInside = CHRINTERVALDICT[chrID];!isInside {
			continue bedLoop
		}

		if _, isInside = PEAKCHRINTERVALDICT[chrID];!isInside {
			continue bedLoop
		}

		processOneRead(&countDicts, inter, cellID, chrID, useGeneCountCompl)
	}

	return countDicts
}


func processOneRead(
	countDicts * map[string]map[string]int,
	inter interval.IntInterface,
	cellID, chrID string,
	useGeneCountCompl bool) {

	var intervalsSNV, intervalsPeaks []interval.IntInterface
	var pos, pos2 int
	var gene, snvID, peakID, eQTLID string
	var startSnv int
	var buffer bytes.Buffer
	var isInside bool

	intervalsSNV = CHRINTERVALDICT[chrID].Get(inter)
	intervalsPeaks = PEAKCHRINTERVALDICT[chrID].Get(inter)

	if len(intervalsSNV) == 0 {
		return
	}

	if len(intervalsPeaks) == 0 {
		return
	}

	for pos = range intervalsSNV {
		gene = GENEUINTDICT[intervalsSNV[pos].ID()]
		startSnv = intervalsSNV[pos].Range().Start
		buffer.WriteString(chrID)
		buffer.WriteRune('_')
		buffer.WriteString(strconv.Itoa(startSnv))
		snvID = buffer.String()
		buffer.Reset()

		if _, isInside = CELLIDEQTLDICT[cellID][snvID];!isInside {
			CELLIDEQTLDICT[cellID][snvID] = true

			for pos2 = range intervalsPeaks {

				peakID = PEAKUNINTDICT[intervalsPeaks[pos2].ID()]

				buffer.WriteString(peakID)
				buffer.WriteRune('_')
				buffer.WriteRune('_')
				buffer.WriteString(gene)

				eQTLID = buffer.String()
				buffer.Reset()

				if useGeneCountCompl {
					(*countDicts)["Score compl."][eQTLID]++
				} else {
					(*countDicts)["Score cluster"][eQTLID]++
				}


				if _, isInside = CELLIDEQTLDICT[cellID][eQTLID];!isInside {
					CELLIDEQTLDICT[cellID][eQTLID] = true

					if useGeneCountCompl {
						(*countDicts)["Cell number compl."][eQTLID]++
					} else {
						(*countDicts)["Cell number cluster"][eQTLID]++
					}
				}
			}
		}
	}
}

func createGeneArray() {
	var line string
	var split []string
	var score float64
	var err error
	var count int

	scanner, file := utils.ReturnReader(GENEFILE, 0)
	defer utils.CloseFile(file)

	GENEDICT = make(map[string]float64)

	for scanner.Scan() {
		line = scanner.Text()
		split = strings.Split(line, "\t")

		if len(split) < 2 {
			log.Fatal(fmt.Sprintf(
				"Error cannot split line: %s in more than 2 part \n", line))
		}

		score, err = strconv.ParseFloat(strings.Trim(split[len(split) - 1], " \n\t\r"), 64)

		if err != nil {
			fmt.Printf("Error with line: %s number: %s\n", line, split[1])
			score = 0

		}

		count++

		GENEDICT[split[0]] = score
	}

	fmt.Printf("Number of genes detected: %d\n", count)

}

func createRefGeneDict() {
	var split []string
	var line string

	if len(GENEDICT) == 0 {
		panic(fmt.Sprintf("Error: GENEDICT must contain at least one gene from %s", GENEFILE))
	}

	EGENEDICT = make(map[string]string)

	if GENEIDTONAMEFILE == "" {
		for gene := range GENEDICT {
			EGENEDICT[gene] = gene
		}
		return
	}

	scanner, file := utils.ReturnReader(GENEIDTONAMEFILE, 0)
	defer utils.CloseFile(file)

	scanner.Scan()

	for scanner.Scan() {
		line = scanner.Text()
		split = strings.Split(line, "\t")

		if _, isInside := GENEDICT[split[1]];!isInside {
			continue
		}

		EGENEDICT[split[0]] = split[1]
	}

}

func createPeakIntervalTree() {
	var line, chrID string
	var split []string
	var err error
	var inter IntInterval
	var peakID bytes.Buffer
	var isInside bool

	var count, start, end int

	PEAKCHRINTERVALDICT = make(map[string]*interval.IntTree)
	PEAKUNINTDICT = make(map[uintptr]string)

	scanner, file := utils.ReturnReader(PEAKFILE, 0)
	defer utils.CloseFile(file)

	for scanner.Scan() {
		line = scanner.Text()
		split = strings.Split(line, "\t")

		chrID = split[0][3:]

		if _, isInside = PEAKCHRINTERVALDICT[chrID];!isInside {
			PEAKCHRINTERVALDICT[chrID] = &interval.IntTree{}
		}

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		inter = IntInterval{
			Start:start,
			End: end,
			UID: uintptr(count),
		}

		peakID.WriteString(split[0])
		peakID.WriteRune('_')
		peakID.WriteString(split[1])
		peakID.WriteRune('_')
		peakID.WriteString(split[2])

		PEAKUNINTDICT[uintptr(count)] = peakID.String()
		peakID.Reset()

		err = PEAKCHRINTERVALDICT[chrID].Insert(inter, false)
		utils.Check(err)

		count++
	}

	fmt.Printf("Number of peak detected: %d\n", count)

}

func createSNPIntervalTree() {
	var line, gene, chrID, fname string
	var split, snv []string
	var err error
	var count, count2, start, end int
	var isInside bool
	var inter IntInterval
	var intervals []interval.IntInterface
	var snpBedFile io.WriteCloser
	var buffer bytes.Buffer

	if WRITESNPTOBED {
		ext := path.Ext(OUTFILE)
		fname = fmt.Sprintf("%s.SNP.bed", OUTFILE[:len(OUTFILE)-len(ext)])
		snpBedFile = utils.ReturnWriter(fname)
		defer utils.CloseFile(snpBedFile)
		defer fmt.Printf("File: %s written\n", fname)
	}

	if len(PEAKCHRINTERVALDICT) == 0 {
		panic(fmt.Sprintf("Error peak intervals seem empty from file: %s\n", PEAKFILE))

	}

	CHRINTERVALDICT = make(map[string]*interval.IntTree)

	geneUintDict := make(map[string]uintptr)
	GENEUINTDICT = make(map[uintptr]string)

	for gene := range GENEDICT {
		geneUintDict[gene] = uintptr(count)
		GENEUINTDICT[uintptr(count)] = gene
		count++
	}

	count = 0

	for _, file := range SNPFILES{

		scanner, file := utils.ReturnReader(file, 0)
		scanner.Scan()

		snpFileLoop:
		for scanner.Scan() {
			line = scanner.Text()
			split = strings.Split(line, "\t")

			if len(split) < 12 {
				log.Fatal(fmt.Sprintf(
					"Error cannot split line: %s in more than 12 part using\t\n", split))
			}

			if _, isInside = EGENEDICT[split[1]];!isInside {
				continue snpFileLoop
			}

			if gene, isInside = EGENEDICT[split[1]];!isInside {
				log.Fatal(fmt.Sprintf("Error gene: %s %s not in ref file:%s\n", split[1], gene, GENEIDTONAMEFILE))

			}

			snv = strings.Split(split[0], "_")

			chrID = snv[0]
			start, err = strconv.Atoi(snv[1])
			utils.Check(err)

			end = start + SNPPADDING
			start = start - SNPPADDING

			count2++

			inter = IntInterval{
				Start:start,
				End: end,
				UID: geneUintDict[gene],
			}

			intervals = PEAKCHRINTERVALDICT[chrID].Get(inter)

			isInside = len(intervals) != 0

			if WRITESNPTOBED {
				buffer.WriteString("chr")
				buffer.WriteString(chrID)
				buffer.WriteRune('\t')
				buffer.WriteString(strconv.Itoa(start))
				buffer.WriteRune('\t')
				buffer.WriteString(strconv.Itoa(end))
				buffer.WriteRune('\t')
				buffer.WriteString(gene)
				buffer.WriteRune('\t')
				buffer.WriteString(strconv.FormatBool(isInside))
				buffer.WriteRune('\n')

				_, err = snpBedFile.Write(buffer.Bytes())
				utils.Check(err)
				buffer.Reset()
			}

			if !isInside {
				continue snpFileLoop
			}

			if _, isInside = CHRINTERVALDICT[chrID];!isInside {
				CHRINTERVALDICT[chrID] = &interval.IntTree{}
			}

			err = CHRINTERVALDICT[chrID].Insert(inter, false)
			utils.Check(err)
			count++
		}

		utils.CloseFile(file)
	}

	fmt.Printf("Number of SNPs linked to significant genes detected: %d\n", count2)
	fmt.Printf("Number of SNPs intersecting with peaks detected: %d\n", count)

}
