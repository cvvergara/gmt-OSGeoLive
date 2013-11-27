#!/bin/bash
#	$Id: GMT_cassini.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt pscoast -R7:30/38:30/10:30/41:30r -JC8.75/40/2.5i -Bafg -Lf9.5/38.8/40/60 -Dh -Gspringgreen \
	-Sazure -Wthinnest -Ia/thinner -P --FONT_LABEL=12p > GMT_cassini.ps
