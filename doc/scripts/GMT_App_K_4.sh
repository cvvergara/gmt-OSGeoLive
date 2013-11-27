#!/bin/bash
#	$Id: GMT_App_K_4.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt pscoast -Rk-100/100/-100/100 -JE130.35/-0.2/3.5i -P -Dh -A1 \
	-Gburlywood -Sazure -Wthinnest -N1/thinnest,- -B30mg10m -BWSne -K > GMT_App_K_4.ps
gmt psbasemap -R -J -O -Dk40+c130.35/-0.2+pthicker >> GMT_App_K_4.ps
