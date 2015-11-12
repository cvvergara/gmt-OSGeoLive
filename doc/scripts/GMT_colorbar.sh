#!/bin/bash
#	$Id: GMT_colorbar.sh 15178 2015-11-06 10:45:03Z fwobbe $
#
ps=GMT_colorbar.ps
gmt makecpt -T-30/30/5 -Cpolar -Z > t.cpt
gmt psbasemap -R0/20/0/1 -JM5i -BWSe -P -K -Baf > $ps
gmt psscale -Ct.cpt -R -J -O -K -Baf -Bx+u"\\232" -By+l@~D@~T -DJBC+o0i/0.35i+w4.5i/0.1i+h+e >> $ps
gmt psxy -R -J -O -T >> $ps
