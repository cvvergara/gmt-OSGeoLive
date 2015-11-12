#!/bin/bash
#	$Id: GMT_images.sh 15178 2015-11-06 10:45:03Z fwobbe $
#
ps=GMT_images.ps
gmt psimage "${src:-.}"/nsf1.jpg -R0/2/0/1 -JX5i/1.6i -B0 -P -K -DjML+w1.5i+o0.1i/0i > $ps
gmt psimage "${src:-.}"/soest.eps -R -J -O -K -DjMR+o0.1i/0i+w2i >> $ps
gmt psxy -R -J -O -T >> $ps
