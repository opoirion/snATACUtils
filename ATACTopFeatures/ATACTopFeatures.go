package main


import (
	"flag"
	utils "gitlab.com/Grouumf/ATACdemultiplex/ATACdemultiplexUtils"
)


/*BEDFILENAME bed file name (input) */
var BEDFILENAME utils.Filename

/*PEAKFILE peak file name (input) */
var PEAKFILE utils.Filename

/*CLUSTERFILE cluster file name (input) */
var CLUSTERFILE utils.Filename

/*FILENAMEOUT  output file name output */
var FILENAMEOUT string


func main() {
	flag.Var(&BEDFILENAME, "bed", "name of the bed file")
	flag.Var(&PEAKFILE, "peak", "File containing peaks")
	flag.Var(&CLUSTERFILE, "cluster", "File containing cluster")
	flag.StringVar(&FILENAMEOUT, "out", "", "name the output file(s)")

	flag.Parse()
}
