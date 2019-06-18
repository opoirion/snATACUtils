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

/*SEP separator used to split the SNV-SNV link file */
var SEP string

/*SNPFILES multiple input files for eQTL SNP-genes */
var SNPFILES utils.ArrayFlags

/*GENEIDTONAMEFILE gene ID to name conversion File */
var GENEIDTONAMEFILE string

/*DBSNPFILE db SNP file  */
var DBSNPFILE string

/*WRITESNPTOBED write significant SNP to bed file */
var WRITESNPTOBED bool

/*USESNPASID Use SNP as id rather than peaks */
var USESNPASID bool

/*CHRINTERVALDICT chr ID <-> interval tree */
var CHRINTERVALDICT map[string]*interval.IntTree

/*PEAKCHRINTERVALDICT chr ID <-> interval tree */
var PEAKCHRINTERVALDICT map[string]*interval.IntTree


/*snpID Internal structure used to represent an SNP */
type snpID struct {
	chrID string
	pos int
}


/*item Internal structure used to sort map */
type item struct {
	key string
	value float64
}

/*SNPTOEXCLUDE snpID -> bool */
var SNPTOEXCLUDE map[snpID]bool

/*SNPDICT chrID <string> -> genomic pos <int> -> score <float64> */
var SNPDICT map[string]map[int]float64

/*SNPIDTOSNPDICT dbSNP ID -> pos score  */
var SNPIDTOSNPDICT map[string]snpID

/*SNPTOSNPIDDICT pos score to dbSNP ID */
var SNPTOSNPIDDICT map[snpID]string

/*GENEDICT genename <string> -> score <float> */
var GENEDICT map[string]float64

/*GENEUINTDICT genename <string> -> score <float> */
var GENEUINTDICT map[uintptr]string

/*PEAKUNINTDICT genename <string> -> score <float> */
var PEAKUNINTDICT map[uintptr]string

/*CELLIDEQTLDICT  <string> -> <string> -> bool */
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

/*CELLIDFNAME file name */
var CELLIDFNAME string

/*CELLIDCOMPLFNAME file name */
var CELLIDCOMPLFNAME string

/*LDLINKFNAME file name */
var LDLINKFNAME string

/*LDFILE file name */
var LDFILE string

/*SNPPADDING Padding for SNP inclusion */
var SNPPADDING int

/*PEAKPADDING Padding for peak */
var PEAKPADDING int

/*PERFORMEQTL perform eQTL analysis */
var PERFORMEQTL bool

/*CREATEDLGROUP create DL group analysis */
var CREATEDLGROUP bool

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
	flag.IntVar(&SNPPADDING, "padding_snp", 0, "Padding added before and after SNP location for SNP inclusion")
	flag.IntVar(&PEAKPADDING, "padding_peak", 0, "Padding added before and after peak coordinates")
	flag.StringVar(&BEDFILENAME, "bed", "", "name of the bed file")
	flag.StringVar(&PEAKFILE, "peak_file", "", "Cluster-specific Peak file ")
	flag.StringVar(&DBSNPFILE, "dbSNP", "", "dbSNP file ")
	flag.StringVar(&CELLIDFNAME, "xgi", "", "cluster specific cell barcodes ")
	flag.StringVar(&CELLIDCOMPLFNAME, "xgi_compl", "", "All cell barcodes ")
	flag.StringVar(&GENEFILE, "gene_file", "", "Cluster-specific Gene file ")
	flag.StringVar(&SEP, "sep", " ", "Separator used to split the SNV-SNV link file ")
	flag.StringVar(&LDFILE, "ld_file", "", "LD file. Each line represents a LD group with SNP sorted using their p-values. The first SNP of each line is considered as the top SNP for each group. The rest of the SNPs will be discarded if found.")
	flag.StringVar(&LDLINKFNAME, "ld_pair", "", "LD file. Each line is a SNV-SNV (<dbSNP1><space><dbSNP2>) link representing a linkage disiquilibrium")
	flag.StringVar(&GENEIDTONAMEFILE, "gene_ID_to_name", "", "File used for gene name conversion (gene_id<tab>gene_name) ")
	flag.Var(&SNPFILES, "eQTL_gene_pair", "name of one or multiple eQTL files (-eQTL_gene_pair <fname1> -eQTL_gene_pair <fname2>)")
	flag.StringVar(&OUTFILE, "out", "", "Name of the output file")
	flag.BoolVar(&USESNPASID, "use_snp_id", false, "use SNP ID rather than peaks as ID")
	flag.BoolVar(&PERFORMEQTL, "analysis", false, `perform eQTL analysis
                        USAGE: ATACeQTLUtils -analysis -xgi <fname> -xgi_compl <fname> -bed <fname> -peak_file <fname> -eQTL_gene_pair <fname> (optionnal: -out <string> -gene_ID_to_name <fname> -ld_file -gene_file <fname>)`)

	flag.BoolVar(&CREATEDLGROUP, "create_dl_group", false, `create DL group
                        USAGE: ATACeQTLUtils -create_dl_group  -eQTL_gene_pair <fname> -dbSNP <fname> -dl_pair <string> (-out <string>)`)

	flag.Parse()

	switch {
	case PERFORMEQTL:
		performEQTLAnalsysis()
	case CREATEDLGROUP:
		createDLGroup()
	}
}


func performEQTLAnalsysis() {

	switch {

	case BEDFILENAME == "":
		log.Fatal(fmt.Printf("!!!! Error -bed option must be provided\n"))
	case PEAKFILE == "":
		log.Fatal(fmt.Printf("!!!! Error -peak option must be provided\n"))
	case CELLIDFNAME == "":
		log.Fatal(fmt.Printf("!!!! Error -xgi must be provided\n"))
	case CELLIDCOMPLFNAME == "":
		log.Fatal(fmt.Printf("!!!! Error -xgi_compl must be provided\n"))
	case len(SNPFILES) == 0:
		log.Fatal(fmt.Printf("!!!! Error at least one -eQTL_gene_pair must be provided\n"))
	case OUTFILE == "":
		OUTFILE = fmt.Sprintf("%s.eQTL.tsv", BEDFILENAME)

	}

	tStart := time.Now()

	if LDFILE != "" {
		//Load LD group
		loadRefLDFile()
	}
	//Create and init cell ID dicts
	loadCellIDDicts()

	if GENEFILE != "" {
		//Create Ref Gene Array
		createGeneArray()
	}

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


func loadRefLDFile() {
	var pos int
	var err error
	var line, chrID, key string
	var split, split2 []string
	var snp  snpID

	snpToKeep := make(map[snpID]bool)
	SNPTOEXCLUDE = make(map[snpID]bool)

	scanner, file := utils.ReturnReader(LDFILE, 0)
	defer utils.CloseFile(file)
	tStart := time.Now()

	for scanner.Scan() {
		line = scanner.Text()
		split = strings.Split(line, "\t")

		split2 = strings.Split(split[0], "_")

		if len(split2) != 2 {
			panic(fmt.Sprintf("Error: SNV: %s\n", split[0]))
		}

		chrID = split2[0]
		pos, err = strconv.Atoi(split2[1])
		utils.Check(err)
		snp.chrID = chrID
		snp.pos = pos

		snpToKeep[snp] = true

		for _, key = range split[1:] {

			split2 = strings.Split(key, "_")

			if len(split2) != 2 {
				panic(fmt.Sprintf("Error: SNV: %s\n", key))
			}

			chrID = split2[0]
			pos, err = strconv.Atoi(split2[1])
			utils.Check(err)
			snp.chrID = chrID
			snp.pos = pos

			SNPTOEXCLUDE[snp] = true
		}
	}

	for snp = range snpToKeep {
		delete(SNPTOEXCLUDE, snp)
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Number of top SNP linked to LD groups %d. Number of excluded SNPs %d \n",
		len(snpToKeep), len(SNPTOEXCLUDE))
	fmt.Printf("Ref DL file loaded in time: %f s \n", tDiff.Seconds())
}


func scanDBSNPFile() {
	var line, chrID, dbSNPID string
	var pos, count int
	var split []string
	var err error
	var isInside bool
	var snp snpID

	scanner, file := utils.ReturnReader(DBSNPFILE, 0)
	defer utils.CloseFile(file)

	tStart := time.Now()

	for scanner.Scan() {
		line = scanner.Text()
		split = strings.Split(line, "\t")

		chrID = split[1][3:]
		pos, err = strconv.Atoi(split[3])

		if err != nil {
			panic(fmt.Sprintf("Error with line: %s and int conversion:%s \n",
				line, split[3]))
		}

		if _,isInside = SNPDICT[chrID]; !isInside {
			continue
		}

		if _,isInside = SNPDICT[chrID][pos]; !isInside {
			continue
		}

		dbSNPID = split[4]
		snp.pos = pos
		snp.chrID = chrID

		SNPIDTOSNPDICT[dbSNPID] = snp
		SNPTOSNPIDDICT[snp] = dbSNPID
		count++
	}

	tDiff := time.Since(tStart)
	fmt.Printf("Scanning dbSNP file done in time: %f s \n", tDiff.Seconds())
	fmt.Printf("Number of sbSNP entries found: %d\n", count)
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
	var startSnv, endSnv int
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

			intervalLoop:
			for pos2 = range intervalsPeaks {

				switch {
				case USESNPASID:
					endSnv = intervalsSNV[pos].Range().End

					buffer.WriteString("chr")
					buffer.WriteString(snvID)
					buffer.WriteRune('_')
					buffer.WriteString(strconv.Itoa(endSnv))
					buffer.WriteString("__")
					buffer.WriteString(gene)

				default:
					peakID = PEAKUNINTDICT[intervalsPeaks[pos2].ID()]

					buffer.WriteString(peakID)
					buffer.WriteRune('_')
					buffer.WriteRune('_')
					buffer.WriteString(gene)
				}

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

				if USESNPASID {
					break intervalLoop
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

	EGENEDICT = make(map[string]string)

	if GENEIDTONAMEFILE == "" {
		for gene := range GENEDICT {
			EGENEDICT[gene] = gene
		}
		return
	}

	isGeneList := GENEFILE != ""
	scanner, file := utils.ReturnReader(GENEIDTONAMEFILE, 0)
	defer utils.CloseFile(file)

	scanner.Scan()

	for scanner.Scan() {
		line = scanner.Text()
		split = strings.Split(line, "\t")

		_, isInside := GENEDICT[split[1]]

		if  isGeneList && !isInside {
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
			Start:start + PEAKPADDING,
			End: end - PEAKPADDING,
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
	var count, count2, count3, start, end int
	var isInside bool
	var inter IntInterval
	var intervals []interval.IntInterface
	var snpBedFile io.WriteCloser
	var buffer bytes.Buffer
	var snp snpID

	checkRefSNP := LDFILE != ""

	isRefGeneList := len(GENEDICT) != 0

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

			gene, isInside = EGENEDICT[split[1]]

			switch {

			case !isInside:
				continue snpFileLoop

			case !isRefGeneList:
				if _, isInside = geneUintDict[gene]; !isInside {
					geneUintDict[gene] = uintptr(count3)
					GENEUINTDICT[uintptr(count3)] = gene
					count3++
				}

			}

			snv = strings.Split(split[0], "_")

			chrID = snv[0]
			start, err = strconv.Atoi(snv[1])
			utils.Check(err)

			if checkRefSNP {
				snp.pos = start
				snp.chrID = chrID

				if SNPTOEXCLUDE[snp] {
					continue snpFileLoop
				}
			}

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
