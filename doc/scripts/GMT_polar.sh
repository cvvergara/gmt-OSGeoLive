#!/bin/bash
#	$Id: GMT_polar.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt grdmath -R0/360/2/4 -I6/0.1 X 4 MUL PI MUL 180 DIV COS Y 2 POW MUL = tt.nc
gmt grdcontour tt.nc -JP3i -B30 -BNs+ghoneydew -P -C2 -S4 --FORMAT_GEO_MAP=+ddd > GMT_polar.ps
