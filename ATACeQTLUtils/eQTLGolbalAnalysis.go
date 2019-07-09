package main


import (
	"log"
	"fmt"
	"strings"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
	"strconv"
	"github.com/biogo/store/interval"
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
	writeGlobaleQTLOutput()
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

func scanPtable() {
	var split []string
	var chro, cluster, gene string
	var start, end int
	var err error
	var intervals []interval.IntInterface
	var interval interval.IntInterface
	var isInside bool

	scanner, file := FEATUREPVALUEFILE.ReturnReader(0)
	defer utils.CloseFile(file)

	//Skip first line
	scanner.Scan()

	for scanner.Scan() {
		split = strings.Split(scanner.Text(), "\t")
		chro = split[0][3:]

		start, err = strconv.Atoi(split[1])
		utils.Check(err)

		end, err = strconv.Atoi(split[2])
		utils.Check(err)

		cluster = split[3]

		if _, isInside = CHRINTERVALDICT[chro];!isInside {
			continue
		}

		intervals = CHRINTERVALDICT[chro].Get(IntInterval{
			Start:start, End:end})

		PEAKRESTCOUNT[cluster]++
		CLUSTERSPECIFICCOUNT[cluster] += len(intervals)

		for _, interval = range intervals {
			gene = GENEUINTDICT[interval.ID()]

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

		if _, isInside = CHRINTERVALDICT[chro];!isInside {
			continue
		}

		intervals = CHRINTERVALDICT[chro].Get(IntInterval{
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
