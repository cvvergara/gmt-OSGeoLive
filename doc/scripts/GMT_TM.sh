#!/bin/bash
#	$Id: GMT_TM.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt pscoast -R0/360/-80/80 -JT330/-45/3.5i -Ba30g -BWSne -Dc -A2000 -Slightblue -G0 -P > GMT_TM.ps
