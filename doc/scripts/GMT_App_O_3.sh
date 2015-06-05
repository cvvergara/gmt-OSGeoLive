#!/bin/bash
#	$Id: GMT_App_O_3.sh 13477 2014-08-27 23:22:00Z pwessel $
#
#	Makes Fig 3 for Appendix O (labeled lines)
#
cat << EOF > fix.txt
80      -8.5
55      -7.5
102     0
130     10.5
EOF
gmt pscoast -R50/160/-15/15 -JM5.3i -Gburlywood -Sazure -A500 -K -P > GMT_App_O_3.ps
gmt grdcontour geoid.nc -J -O -B20f10 -BWSne -C10 -A20+d+f8p -Gffix.txt/0.1i -S10 -T:LH >> GMT_App_O_3.ps
