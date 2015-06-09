#!/bin/bash
#
#	$Id: annotations1.sh 12114 2013-09-03 19:19:00Z fwobbe $

ps=annotations1.ps

basemap="gmt psbasemap -JX3id/2.5id --FONT_ANNOT_PRIMARY=10p"
$basemap -R-25/25/-19/23 -B10 -BWSne -P -K --FORMAT_GEO_MAP=dddF -Xf1i -Yf1i > $ps
$basemap -R-1.5/1.5/-1.2/1.5 -B0.5 -BwSnE -O -K  --FORMAT_GEO_MAP=ddd.xF -Xf4.5i >> $ps
$basemap -R-1.05/1.05/-1.1/1.3 -B30m -BWSne -O -K --FORMAT_GEO_MAP=ddd:mmF -Xf1i -Yf4.25i >> $ps
$basemap -R-0:00:50/0:01:00/-0:01:00/0:01:00 -B0.5m -BwSnE -O -K  --FORMAT_GEO_MAP=ddd:mm.xF -Xf4.5i >> $ps
$basemap -R-0:00:30/0:00:30/-0:01:00/0:01:00 -B30s -BWSne -O -K --FORMAT_GEO_MAP=ddd:mm:ssF -Xf1i -Yf7.5i >> $ps
$basemap -R-0:00:04/0:00:05/-0:00:05/0:00:05 -B2.5s -BwSnE -O --FORMAT_GEO_MAP=ddd:mm:ss.xF -Xf4.5i >> $ps
