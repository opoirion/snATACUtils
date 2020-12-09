module github.com/opoirion/snATACUtils/BAMutils

replace github.com/opoirion/snATACUtils/ATACdemultiplexUtils => ../ATACdemultiplexUtils

go 1.15

require (
	github.com/biogo/hts v1.2.2
	github.com/opoirion/snATACUtils/ATACdemultiplexUtils v0.0.0-00010101000000-000000000000
)
