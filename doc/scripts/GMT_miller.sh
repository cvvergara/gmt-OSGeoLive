#!/bin/bash
#	$Id: GMT_miller.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt pscoast -R-90/270/-80/90 -Jj1:400000000 -Bx45g45 -By30g30 -Dc -A10000 -Gkhaki -Wthinnest -P \
	-Sazure > GMT_miller.ps
