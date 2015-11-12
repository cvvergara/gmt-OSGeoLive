#!/bin/bash
#	$Id: GMT_tut_4.sh 14523 2015-07-15 02:26:00Z pwessel $
#
gmt pscoast -R-130/-70/24/52 -JB-100/35/33/45/6i -Ba -B+t"Conic Projection" -N1/thickest -N2/thinnest -A500 -Ggray -Wthinnest -P > GMT_tut_4.ps
