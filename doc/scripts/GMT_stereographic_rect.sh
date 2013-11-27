#!/bin/bash
#	$Id: GMT_stereographic_rect.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt gmtset MAP_ANNOT_OBLIQUE 30
gmt pscoast -R-25/59/70/72r -JS10/90/11c -B20g -Dl -A250 -Gdarkbrown -Wthinnest -P \
	-Slightgray > GMT_stereographic_rect.ps
