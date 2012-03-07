#!/bin/csh
#		$Id: job23.csh,v 1.22 2011/03/01 01:34:48 remko Exp $
#		GMT EXAMPLE 22
#
# Purpose:	Plot distances from Rome and draw shortest paths
# GMT progs:	grdmath, grdcontour, pscoast, psxy, pstext, grdtrack
# Unix progs:	echo, cat, awk

# Position and name of central point:

set lon = 12.50
set lat = 41.99
set name = "Rome"

# Calculate distances (km) to all points on a global 1x1 grid

grdmath -Rg -I1 $lon $lat SDIST 111.13 MUL = dist.nc

# Location info for 5 other cities + label justification

cat << END >! cities.d
105.87	21.02	HANOI		LM
282.95	-12.1	LIMA		LM
178.42	-18.13	SUVA		LM
237.67	47.58	SEATTLE		RM
28.20	-25.75	PRETORIA	LM
END

pscoast -Rg -JH90/9i -Glightgreen -Sblue -U"Example 23 in Cookbook" -A1000 \
  -B0g30:."Distances from $name to the World": -K -Dc -Wthinnest >! ../example_23.ps

grdcontour dist.nc -A1000+v+ukm+kwhite -Glz-/z+ -S8 -C500 -O -K -J -Wathin,white \
  -Wcthinnest,white,- >> ../example_23.ps

# Find the number of cities:

set n = `cat cities.d | wc -l`

# For each of these cities, plot great circle arc with psxy

set i = 1
while ($i <= $n)
	set record = `$AWK '{ if (NR == '$i') print $0}' cities.d`
	(echo $lon $lat; echo $record[1] $record[2]) | \
		psxy -R -J -O -K -Wthickest/red >> ../example_23.ps
	@ i++
end

# Plot red squares at cities and plot names:

psxy -R -J -O -K -Ss0.2 -Gred -Wthinnest cities.d >> ../example_23.ps
$AWK '{print $1, $2, 12, 0, 9, $4, $3}' cities.d \
	| pstext -R -J -O -K -Dj0.15/0 -Gred -N >> ../example_23.ps

# Place a yellow star at Rome

echo "$lon $lat" | psxy -R -J -O -K -Sa0.2i -Gyellow -Wthin >> ../example_23.ps

# Sample the distance grid at the cities and use the distance in km for labels

grdtrack -Gdist.nc cities.d \
	| $AWK '{printf "%s %s 12 0 1 CT %d\n", $1, $2, int($NF+0.5)}' \
	| pstext -R -J -O -D0/-0.2i -N -Wwhite,o -C0.02i/0.02i >> ../example_23.ps

# Clean up after ourselves:

\rm -f cities.d dist.nc .gmt*
