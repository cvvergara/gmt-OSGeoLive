#!/bin/bash
#	$Id: GMT_stereographic_polar.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt pscoast -R-30/30/60/72 -Js0/90/4.5i/60 -B10g -Dl -A250 -Groyalblue \
	-Sseashell -P > GMT_stereographic_polar.ps
