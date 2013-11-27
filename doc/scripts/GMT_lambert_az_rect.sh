#!/bin/bash
#	$Id: GMT_lambert_az_rect.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt gmtset FORMAT_GEO_MAP ddd:mm:ssF MAP_GRID_CROSS_SIZE_PRIMARY 0
gmt pscoast -R0/-40/60/-10r -JA30/-30/4.5i -Bag -Dl -A500 -Gp300/10 -Wthinnest -P \
	> GMT_lambert_az_rect.ps
