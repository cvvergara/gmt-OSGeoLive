#!/bin/bash
#	$Id: GMT_perspective.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt pscoast -Rg -JG4/52/230/90/60/180/60/60/5i -Bx2g2 -By1g1 -Ia -Di -Glightbrown -Wthinnest -P \
	-Slightblue --MAP_ANNOT_MIN_SPACING=0.25i > GMT_perspective.ps
