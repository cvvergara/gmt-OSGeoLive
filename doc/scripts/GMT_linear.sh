#!/bin/bash
#	$Id: GMT_linear.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt psxy -R0/100/0/10 -JX3i/1.5i -Bag -BWSne+gsnow -Wthick,blue,- -P -K sqrt.d > GMT_linear.ps
gmt psxy -R -J -St0.1i -N -Gred -Wfaint -O sqrt.d10 >> GMT_linear.ps
