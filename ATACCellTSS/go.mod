module github.com/opoirion/snATACUtils/ATACCellTSS

replace github.com/opoirion/snATACUtils/ATACdemultiplexUtils => ../ATACdemultiplexUtils

go 1.15

require (
	github.com/biogo/store v0.0.0-20201120204734-aad293a2328f
	github.com/opoirion/snATACUtils/ATACdemultiplexUtils v0.0.0-00010101000000-000000000000
)
