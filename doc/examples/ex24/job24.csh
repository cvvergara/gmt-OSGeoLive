#!/bin/csh
#		$Id: job24.csh,v 1.12 2011/03/01 01:34:48 remko Exp $
#		GMT EXAMPLE 24
#
# Purpose:	Extract subsets of data based on geospatial criteria
# GMT progs:	gmtselect, pscoast, psxy, minmax
# Unix progs:	echo, cat, awk

# Highlight oceanic earthquakes within 3000 km of Hobart and > 1000 km from dateline

echo "147:13 -42:48 3000 Hobart" >! point.d
cat << END >! dateline.d
> Our proxy for the dateline
180	0
180	-90
END
set R = `minmax -I10 oz_quakes.d`
pscoast $R -JM9i -K -Gtan -Sdarkblue -Wthin,white -Dl -A500 -Ba20f10g10WeSn \
	-U"Example 24 in Cookbook" >! ../example_24.ps
psxy -R -J -O -K oz_quakes.d -Sc0.05i -Gred >> ../example_24.ps
gmtselect oz_quakes.d -L1000/dateline.d -Nk/s -C3000/point.d -fg -R -Il \
	| psxy -R -JM -O -K -Sc0.05i -Ggreen >> ../example_24.ps
$AWK '{print $1, $2, 0, 2*$3, 2*$3}' point.d | psxy -R -J -O -K -SE -Wfat,white >> ../example_24.ps
$AWK '{print $1, $2, 14, 0, 1, "LT", $4}' point.d \
	| pstext -R -J -O -K -Gwhite -D0.1i/-0.1i >> ../example_24.ps
psxy -R -J -O -K point.d -Wfat,white -S+0.2i >> ../example_24.ps
psxy -R -J -O -m dateline.d -Wfat,white -A >> ../example_24.ps

\rm -f point.d dateline.d .gmt*
