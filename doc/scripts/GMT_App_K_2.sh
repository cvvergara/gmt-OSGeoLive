#!/bin/bash
#	$Id: GMT_App_K_2.sh 15178 2015-11-06 10:45:03Z fwobbe $
#
gmt pscoast -Rk-2000/2000/-2000/2000 -JE130.35/-0.2/3.5i -P -Dl -A100 \
	-Gburlywood -Sazure -Wthinnest -N1/thinnest,- -B10g5 -BWSne -K > GMT_App_K_2.ps
gmt psbasemap -R -J -O -Dg130.35/-0.2+w1000k+jCM -F+pthicker >> GMT_App_K_2.ps
