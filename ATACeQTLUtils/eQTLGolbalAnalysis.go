package main


import (
	"log"
	"fmt"
	"strings"
	"github.com/biogo/store/interval"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"strconv"
	"path"
	"bytes"
	stats "github.com/glycerine/golang-fisher-exact"
	"sort"
)



/*CLUSTERGENEFILE file name */
var CLUSTERGENEFILE utils.Filename

/*CLUSTERRNASPECIFICEQTL <cluster> count*/
var CLUSTERRNASPECIFICEQTL map[string]int

/*EQTLRESTCOUNT <cluster> count*/
var EQTLRESTCOUNT map[string]int

/*PEAKRESTOUTSIDECOUNT <cluster> count*/
var PEAKRESTOUTSIDECOUNT map[string]int

/*PEAKRESTCOUNT <cluster> count*/
var PEAKRESTCOUNT map[string]int

/*CLUSTERSPECIFICCOUNT <cluster> count*/
var CLUSTERSPECIFICCOUNT map[string]int

/*TOTALRNAEQTLCOUNT <cluster> count*/
var TOTALRNAEQTLCOUNT map[string]int

/*CLUSTERGENEDICT map[clusterID]map[gene]bool*/
var CLUSTERGENEDICT map[string]map[string]bool

/*TOTALNUMBEREQTL ...*/
var TOTALNUMBEREQTL int

/*EQTLPEAKRNASPECIFIC <cluster> count*/
var EQTLPEAKRNASPECIFIC  map[string]int

/*EQTLPEAKSPECIFIC <cluster> count*/
var EQTLPEAKSPECIFIC  map[string]int

/*EQTLRNASPECIFIC <cluster> count*/
var EQTLRNASPECIFIC  map[string]int

/*EQTLREST <cluster> count*/
var EQTLREST  map[string]int


func performGlobalEQTLAnalsysis() {

	switch {
	case CLUSTERGENEFILE == "":
		log.Fatal(fmt.Printf("!!!! Error -cluster_gene option must be provided\n"))

	case GENEIDTONAMEFILE == "":
		log.Fatal(fmt.Printf("!!!! Error -gene_ID_to_name option must be provided\n"))

	case FEATUREPVALUEFILE == "":
		log.Fatal(fmt.Printf("!!!! Error -ptable option must be provided\n"))

	case len(SNPFILES) == 0:
		log.Fatal(fmt.Printf("!!!! Error -eQTL option must be provided\n"))

	case PEAKFILE == "":
		log.Fatal(fmt.Printf("!!!! Error -peak option must be provided\n"))

	case OUTFILE == "":
		ext := path.Ext(FEATUREPVALUEFILE.String())
		OUTFILE = fmt.Sprintf("%s.globaleQTL.tsv",
			FEATUREPVALUEFILE[:len(FEATUREPVALUEFILE) - len(ext)])
	}

	PEAKRESTOUTSIDECOUNT = make(map[string]int)
	PEAKRESTCOUNT = make(map[string]int)
	EQTLRESTCOUNT = make(map[string]int)
	CLUSTERRNASPECIFICEQTL = make(map[string]int)
	CLUSTERSPECIFICCOUNT = make(map[string]int)
	TOTALRNAEQTLCOUNT = make(map[string]int)

	if LDFILE != "" {
		//Load LD group
		loadRefLDFile()
	}

	createRefGeneDict()
	loadClusterGeneDict()
	createPeakIntervalTree()
	createSNPIntervalTree()
	createRefGeneDict()
	scanPtable()
	scanPeaks()
	countForNonRNAeQTLCount()
	scanSNPList()
	writeGlobaleQTLOutput()
	writeGlobaleQTLOutputTable2()
}


func countForNonRNAeQTLCount() {
	var gene, cluster string

	for gene = range GENETOEQTLCOUNT {
		TOTALNUMBEREQTL += GENETOEQTLCOUNT[gene]
	}

	for cluster = range CLUSTERGENEDICT {
		for gene = range CLUSTERGENEDICT[cluster] {
			TOTALRNAEQTLCOUNT[cluster] += GENETOEQTLCOUNT[gene]
		}
	}
}

func writeGlobaleQTLOutput() {
	writer := utils.ReturnWriter(OUTFILE)
	defer utils.CloseFile(writer)
	var buffer bytes.Buffer
	var n11, n12, n21, n22 int

	buffer.WriteString(
		"ClusterID\tCluster-RNA-eQTL\tCluster-eQTL\tpeak-EQTL\tnon-cluster-eQTL\tchi2 p-value\n")

	clusterList := []string{}

	for cluster := range CLUSTERGENEDICT {
		clusterList = append(clusterList, strings.Trim(cluster, " "))
	}

	sort.Strings(clusterList)

	for _, cluster := range clusterList {
		n11 = CLUSTERRNASPECIFICEQTL[cluster]
		n12 = TOTALRNAEQTLCOUNT[cluster]
		n21 = CLUSTERSPECIFICCOUNT[cluster] - n11
		n22 = TOTALNUMBEREQTL - n12

		_, pvalue :=  stats.ChiSquareTest(n11, n12,
			n21, n22, true)

		buffer.WriteString(fmt.Sprintf("%s\t%d\t%d\t%d\t%d\t%f\n",
			cluster, n11, n12, n21, n22, pvalue))

	}

	_, err := writer.Write(buffer.Bytes())
	utils.Check(err)

	fmt.Printf("results:\n%s\n", buffer.String())
}


func writeGlobaleQTLOutputTable2() {
	ext := path.Ext(FEATUREPVALUEFILE.String())
	outfile := fmt.Sprintf("%s.globaleQTL.Table2.tsv",
		FEATUREPVALUEFILE[:len(FEATUREPVALUEFILE) - len(ext)])


	writer := utils.ReturnWriter(outfile)
	defer utils.CloseFile(writer)
	var buffer bytes.Buffer
	var n11, n12, n21, n22 int


	fmt.Printf("\n\n####Analysis 2 ####\n")

	buffer.WriteString(
		"ClusterID\tC-specific scRNA + chromatin\tC-specific chromatin\tC-specific scRNA\tOther eQTL\tP-value\n")

	clusterList := []string{}

	for cluster := range CLUSTERGENEDICT {
		clusterList = append(clusterList, strings.Trim(cluster, " "))
	}

	sort.Strings(clusterList)

	for _, cluster := range clusterList {
		n11 = EQTLPEAKRNASPECIFIC[cluster]
		n12 = EQTLPEAKSPECIFIC[cluster]
		n21 = EQTLRNASPECIFIC[cluster]
		n22 = EQTLREST[cluster]

		_, pvalue :=  stats.ChiSquareTest(n11, n12,
			n21, n22, true)

		buffer.WriteString(fmt.Sprintf("%s\t%d\t%d\t%d\t%d\t%f\n",
			cluster, n11, n12, n21, n22, pvalue))

	}

	_, err := writer.Write(buffer.Bytes())
	utils.Check(err)

	fmt.Printf("results:\n%s\n", buffer.String())

	fmt.Printf("\nFile: %s written\n", outfile)
}


func scanPtable() {
	var split []string
	var chro, cluster, gene string
	var start, end int
	var err error
	var intervals []interval.IntInterface
	var inter interval.IntInterface
	var isInside bool
	var snp Snp
	var interRange interval.IntRange

	scanner, file := FEATUREPVALUEFILE.ReturnReader(0)
	defer utils.CloseFile(file)

	//Skip first line
	scanner.Scan()

	SNPCLUSTERLIST = make(map[string]map[Snp]bool)

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), "\t")
		chro = split[0][3:]

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		cluster = split[3]

		if _, isInside = SNPINTERVALDICT[chro];!isInside {
			continue
		}

		if _, isInside = SNPCLUSTERLIST[cluster];!isInside {
			SNPCLUSTERLIST[cluster] = make(map[Snp]bool)
		}

		intervals = SNPINTERVALDICT[chro].Get(IntInterval{
			Start:start, End:end})

		PEAKRESTCOUNT[cluster]++
		CLUSTERSPECIFICCOUNT[cluster] += len(intervals)

		for _, inter = range intervals {
			gene = GENEUINTDICT[inter.ID()]

			interRange = inter.Range()

			snp.chrID = chro
			snp.Start = interRange.Start
			snp.End = interRange.End
			snp.gene = gene

			SNPCLUSTERLIST[cluster][snp] = true

			if CLUSTERGENEDICT[cluster][gene] {
				CLUSTERRNASPECIFICEQTL[cluster]++
			}
		}
	}
}


func scanPeaks() {
	var split []string
	var start, end int
	var chro, cluster, gene string

	var intervals []interval.IntInterface
	var interval interval.IntInterface
	var isInside bool
	var err error

	scanner, file := PEAKFILE.ReturnReader(0)
	defer utils.CloseFile(file)

	for scanner.Scan() {

		split = strings.Split(scanner.Text(), "\t")
		chro = split[0][3:]

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		if _, isInside = SNPINTERVALDICT[chro];!isInside {
			continue
		}

		intervals = SNPINTERVALDICT[chro].Get(IntInterval{
			Start:start, End:end})

		for cluster = range CLUSTERGENEDICT {
			cluster = strings.Trim(cluster, " ")
			PEAKRESTOUTSIDECOUNT[cluster]++

			for _, interval = range intervals {
				gene = GENEUINTDICT[interval.ID()]

				switch cluster {
				case cluster:
					if CLUSTERGENEDICT[cluster][gene] {
						EQTLRESTCOUNT[cluster]++
					}

				}

			}

		}

	}

}


func scanSNPList() {

	EQTLPEAKRNASPECIFIC = make(map[string]int)
	EQTLPEAKSPECIFIC = make(map[string]int)
	EQTLRNASPECIFIC = make(map[string]int)
	EQTLREST = make(map[string]int)

	var snp Snp
	var cluster string

	for snp = range SNPLIST {

		for cluster = range CLUSTERGENEDICT {
			if SNPCLUSTERLIST[cluster][snp] {
				if CLUSTERGENEDICT[cluster][snp.gene] {
					EQTLPEAKRNASPECIFIC[cluster]++
				} else {
					EQTLPEAKSPECIFIC[cluster]++

				}
			} else {
				if CLUSTERGENEDICT[cluster][snp.gene] {
					EQTLRNASPECIFIC[cluster]++
				} else {
					EQTLREST[cluster]++
				}
			}
		}
	}
}


func loadClusterGeneDict() {
	var line, cluster string
	var fname utils.Filename
	var split []string

	CLUSTERGENEDICT = make(map[string]map[string]bool)

	scanner, file := CLUSTERGENEFILE.ReturnReader(0)
	defer utils.CloseFile(file)

	for scanner.Scan() {
		line = scanner.Text()
		split = strings.Split(line, "\t")

		if len(split) < 2 {
			log.Fatal(fmt.Sprintf("error with line %s from file %s. Cannot split with <tab>!\n",
				line, CLUSTERGENEFILE))
		}

		cluster, fname = split[0], utils.Filename(split[1])
		cluster = strings.Trim(cluster, " ")

		loadOneGeneFile(fname, cluster)
	}
}


func loadOneGeneFile(fname utils.Filename, cluster string) {
	var gene string
	CLUSTERGENEDICT[cluster] = make(map[string]bool)

	scanner, file := fname.ReturnReader(0)
	defer utils.CloseFile(file)

	for scanner.Scan() {
		gene = strings.Split(scanner.Text(), "\t")[0]

		CLUSTERGENEDICT[cluster][gene] = true
	}
}
