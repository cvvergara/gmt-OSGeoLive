#!/bin/bash
#	$Id: GMT_App_O_1.sh 11490 2013-05-16 06:26:21Z pwessel $
#
#	Makes Fig 1 for Appendix O (labeled lines)
#
gmt pscoast -R50/160/-15/15 -JM5.3i -Gburlywood -Sazure -A500 -K -P > GMT_App_O_1.ps
gmt grdcontour geoid.nc -J -O -B20f10 -BWSne -C10 -A20+f8p -Gd1.5i -S10 -T:LH >> GMT_App_O_1.ps
