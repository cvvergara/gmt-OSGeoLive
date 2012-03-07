#!/bin/csh
#		$Id: job08.csh,v 1.13 2011/03/01 01:34:48 remko Exp $
#		GMT EXAMPLE 08
#
# Purpose:	Make a 3-D bar plot
# GMT progs:	grd2xyz, pstext, psxyz
# Unix progs:	echo, rm

grd2xyz guinea_bay.nc | \
	psxyz -B1/1/1000:"Topography (m)"::.ETOPO5:WSneZ+ -R-0.1/5.1/-0.1/5.1/-5000/0 \
	-P -JM5i -JZ6i -E200/30 -So0.0833333ub-5000 -U"Example 8 in Cookbook" \
	-Wthinnest -Glightgray -K >! ../example_08.ps
echo '0.1 4.9 24 0 1 TL This is the surface of cube' | pstext -R -J -JZ -Z0 -E200/30 -O \
	>> ../example_08.ps

\rm -f .gmt*
