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

/*EQTLCOUNT <eQTL> count*/
var EQTLCOUNT map[string]int

/*EQTLRESTCOUNT <eQTL> count*/
var EQTLRESTCOUNT map[string]int

/*PEAKRESTOUTSIDECOUNT <eQTL> count*/
var PEAKRESTOUTSIDECOUNT map[string]int

/*PEAKRESTCOUNT <eQTL> count*/
var PEAKRESTCOUNT map[string]int

/*CLUSTERGENEDICT map[clusterID]map[gene]bool*/
var CLUSTERGENEDICT map[string]map[string]bool




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
	EQTLCOUNT = make(map[string]int)

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
	writeGlobaleQTLOutput()
}


func writeGlobaleQTLOutput() {
	writer := utils.ReturnWriter(OUTFILE)
	defer utils.CloseFile(writer)
	var buffer bytes.Buffer

	buffer.WriteString(
		"ClusterID\tCluster-eQTL\tCluster-Peak\tBackground-EQTL\tBackground-Peak\tchi2 p-value\n")

	clusterList := []string{}

	for cluster := range CLUSTERGENEDICT {
		clusterList = append(clusterList, strings.Trim(cluster, " "))
	}

	sort.Strings(clusterList)

	for _, cluster := range clusterList {
		_, pvalue :=  stats.ChiSquareTest(EQTLCOUNT[cluster], PEAKRESTCOUNT[cluster],
			EQTLRESTCOUNT[cluster], PEAKRESTOUTSIDECOUNT[cluster], true)

		buffer.WriteString(fmt.Sprintf("%s\t%d\t%d\t%d\t%d\t%f\n",
			cluster, EQTLCOUNT[cluster], PEAKRESTCOUNT[cluster],
			EQTLRESTCOUNT[cluster], PEAKRESTOUTSIDECOUNT[cluster], pvalue))

	}

	_, err := writer.Write(buffer.Bytes())
	utils.Check(err)

	fmt.Printf("results:\n%s\n", buffer.String())
}

func scanPtable() {
	var split []string
	var chro, cluster, cluster2, gene string
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

		for cluster2 = range CLUSTERGENEDICT {
			cluster2 = strings.Trim(cluster2, " ")
			for _, interval = range intervals {
				gene = GENEUINTDICT[interval.ID()]

				switch cluster2 {
				case cluster:
					if CLUSTERGENEDICT[cluster][gene] {
						EQTLCOUNT[cluster]++
					}

					PEAKRESTCOUNT[cluster]++
				default:
					if CLUSTERGENEDICT[cluster][gene] {
						EQTLRESTCOUNT[cluster]++
					}

					PEAKRESTOUTSIDECOUNT[cluster]++
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
