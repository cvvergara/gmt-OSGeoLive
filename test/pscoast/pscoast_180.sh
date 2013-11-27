#!/bin/bash
#	$Id: pscoast_180.sh 11490 2013-05-16 06:26:21Z pwessel $
#

ps=pscoast_180.ps

gmt pscoast -R-10/180/20/50 -JM6i -Bxa20g20 -Bya20g20+10 -Dc -G221/204/170 -A0/0/1 -P -Xc -K > $ps
gmt pscoast -R0/180/20/50 -JM6i -Bxa20g20 -Bya20g20+10 -Dc -G221/204/170 -A0/0/1 -O -Y2.5i >> $ps
