#!/bin/bash
#	$Id: GMT_log.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt psxy -R1/100/0/10 -Jx1.5il/0.15i -Bx2g3 -Bya2f1g2 -BWSne+gbisque -Wthick,blue,- -P -K -h sqrt.d \
	> GMT_log.ps
gmt psxy -R -J -Ss0.1i -N -Gred -W -O -h sqrt.d10 >> GMT_log.ps
