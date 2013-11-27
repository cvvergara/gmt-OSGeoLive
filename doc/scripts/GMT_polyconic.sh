#!/bin/bash
#	$Id: GMT_polyconic.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt pscoast -R-180/-20/0/90 -JPoly/4i -Bx30g10 -By10g10 -Dc -A1000 -Glightgray -Wthinnest -P \
	> GMT_polyconic.ps
