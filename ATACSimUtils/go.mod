module github.com/opoirion/snATACUtils/ATACSimUtils

replace github.com/opoirion/snATACUtils/ATACdemultiplexUtils => ../ATACdemultiplexUtils

go 1.15

require (
	github.com/valyala/fastrand v1.0.0
	github.com/opoirion/snATACUtils/ATACdemultiplexUtils v0.0.0-00010101000000-000000000000
)
