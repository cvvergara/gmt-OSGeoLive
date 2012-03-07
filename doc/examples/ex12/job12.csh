#!/bin/csh
#		$Id: job12.csh,v 1.14 2011/03/01 01:34:48 remko Exp $
#		GMT EXAMPLE 12
#
# Purpose:	Illustrates Delaunay triangulation of points, and contouring
# GMT progs:	makecpt, minmax, pscontour, pstext, psxy, triangulate
# Unix progs:	$AWK, echo, rm

# First draw network and label the nodes

triangulate table_5.11 -m >! net.xy
psxy -R0/6.5/-0.2/6.5 -JX3.06i/3.15i -B2f1WSNe -m net.xy -Wthinner -P -K -X0.9i -Y4.65i \
	>! ../example_12.ps
psxy table_5.11 -R -J -O -K -Sc0.12i -Gwhite -Wthinnest >> ../example_12.ps
$AWK '{print $1, $2, 6, 0, 0, "CM", NR-1}' table_5.11 | \
	pstext -R -J -O -K >> ../example_12.ps

# Then draw network and print the node values

psxy -R -J -B2f1eSNw -m net.xy -Wthinner -O -K -X3.25i >> ../example_12.ps
psxy -R -J -O -K table_5.11 -Sc0.03i -Gblack >> ../example_12.ps
$AWK '{printf "%g %s 6 0 0 LM %g\n", $1, $2, $3}' table_5.11 | pstext -R -J -O -K -Wwhite,o \
	-C0.01i/0.01i -D0.08i/0i -N >> ../example_12.ps

# Then contour the data and draw triangles using dashed pen
# Use "minmax" and "makecpt" to make a color palette (.cpt) file

set T = `minmax -T25/2 table_5.11`
makecpt -Cjet $T >! topo.cpt
pscontour -R -J table_5.11 -B2f1WSne -Wthin -Ctopo.cpt -Lthinnest,- -G1i/0 -X-3.25i -Y-3.65i -O -K \
	-U"Example 12 in Cookbook" >> ../example_12.ps

# Finally color the topography

pscontour -R -J table_5.11 -B2f1eSnw -Ctopo.cpt -I -X3.25i -O -K >> ../example_12.ps
echo "3.16 8 30 0 1 BC Delaunay Triangulation" | pstext -R0/8/0/11 -Jx1i -O -X-3.25i >> ../example_12.ps

\rm -f net.xy topo.cpt .gmt*
