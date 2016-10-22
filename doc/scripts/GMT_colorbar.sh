#!/bin/bash
#	$Id: GMT_colorbar.sh 16573 2016-06-19 02:42:45Z pwessel $
#
ps=GMT_colorbar.ps
gmt makecpt -T-15/15 -Cpolar > t.cpt
gmt psbasemap -R0/20/0/1 -JM5i -BWSe -P -K -Baf > $ps
gmt psscale -Ct.cpt -R -J -O -K -Baf -Bx+u"\\232" -By+l@~D@~T -DJBC+o0i/0.35i+w4.5i/0.1i+h+e >> $ps
gmt psxy -R -J -O -T >> $ps
