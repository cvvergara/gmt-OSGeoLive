#!/bin/csh
#		$Id: job15.csh 9545 2011-07-27 19:31:54Z pwessel $
#		GMT EXAMPLE 15
#
# Purpose:	Gridding and clipping when data are missing
# GMT progs:	blockmedian, gmtconvert, grdclip, grdcontour, grdinfo, minmax
# GMT progs:	nearneighbor, pscoast, psmask, pstext, surface
# Unix progs:	echo, rm

gmtconvert ship.xyz -bo >! ship.b
set region = `minmax ship.b -I1 -bi3`
nearneighbor $region -I10m -S40k -Gship.nc ship.b -bi3
set info = `grdinfo -C -M ship.nc`
grdcontour ship.nc -JM3i -P -B2WSne -C250 -A1000 -G2i -K -U"Example 15 in Cookbook" >! ../example_15.ps

blockmedian $region -I10m ship.b -bi3 -bo >! ship_10m.b
surface $region -I10m ship_10m.b -Gship.nc -bi3
psmask $region -I10m ship.b -J -O -K -T -Glightgray -bi3 -X3.6i >> ../example_15.ps
grdcontour ship.nc -J -B2WSne -C250 -L-8000/0 -A1000 -G2i -O -K >> ../example_15.ps

psmask $region -I10m ship_10m.b -bi3 -J -B2WSne -O -K -X-3.6i -Y3.75i >> ../example_15.ps
grdcontour ship.nc -J -C250 -A1000 -L-8000/0 -G2i -O -K >> ../example_15.ps
psmask -C -O -K >> ../example_15.ps

grdclip ship.nc -Sa-1/NaN -Gship_clipped.nc
grdcontour ship_clipped.nc -J -B2WSne -C250 -A1000 -L-8000/0 -G2i -O -K -X3.6i >> ../example_15.ps
pscoast $region -J -O -K -Ggray -Wthinnest >> ../example_15.ps
echo $info[12] $info[13] | psxy -R -J -O -K -Sa0.15i -Wthick >> ../example_15.ps
echo "-0.3 3.6 24 0 1 CB Gridding with missing data" | pstext -R0/3/0/4 -Jx1i -O -N >> ../example_15.ps

\rm -f ship.b ship_10m.b ship.nc ship_clipped.nc .gmt*
