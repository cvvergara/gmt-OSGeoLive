#!/bin/bash
#	$Id: GMT_gnomonic.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt pscoast -Rg -JF-120/35/60/4.5i -B30g15 -Dc -A10000 -Gtan -Scyan -Wthinnest -P > GMT_gnomonic.ps
