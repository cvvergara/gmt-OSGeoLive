#!/bin/bash
#
#	$Id: annotations4.sh 12114 2013-09-03 19:19:00Z fwobbe $

ps=annotations4.ps

gmt psbasemap -R-180/180/-90/90 -B60g60 -JX14cd/7cd --MAP_FRAME_TYPE=inside --FONT_ANNOT_PRIMARY=10p -P -K > $ps
gmt psbasemap -R-180/180/-90/90 -B60g60 -JX14cd/7cd --MAP_FRAME_TYPE=plain  --FONT_ANNOT_PRIMARY=10p -O >> $ps
