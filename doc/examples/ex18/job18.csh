#!/bin/csh
#		$Id: job18.csh 9606 2011-10-19 17:52:57Z pwessel $
#		GMT EXAMPLE 18
#
# Purpose:	Illustrates volumes of grids inside contours and spatial
#		selection of data
# GMT progs:	gmtset, gmtselect, grdclip, grdcontour, grdgradient, grdimage
# GMT progs:	grdmath, grdvolume, makecpt, pscoast, psscale, pstext, psxy
# Unix progs:	$AWK, cat, rm

# Use spherical projection since SS data define on sphere

gmtset ELLIPSOID Sphere D_FORMAT %g

# Define location of Pratt seamount

echo "-142.65	56.25" >! pratt.d

# First generate gravity image w/ shading, label Pratt, and draw a circle
# of radius = 200 km centered on Pratt.

makecpt -Crainbow -T-60/60/10 -Z >! grav.cpt
grdgradient AK_gulf_grav.nc -Nt1 -A45 -GAK_gulf_grav_i.nc
grdimage AK_gulf_grav.nc -IAK_gulf_grav_i.nc -JM5.5i -Cgrav.cpt -B2f1 -P -K -X1.5i -Y5.85i \
	>! ../example_18.ps
pscoast -RAK_gulf_grav.nc -J -O -K -Di -Ggray -Wthinnest >> ../example_18.ps
psscale -D2.75i/-0.4i/4i/0.15ih -Cgrav.cpt -B20f10/:mGal: -O -K >> ../example_18.ps
$AWK '{print $1, $2, 12, 0, 1, "LB", "Pratt"}' pratt.d | pstext -R -J -O -K -D0.1i/0.1i \
	>> ../example_18.ps
$AWK '{print $1, $2, 0, 400, 400}' pratt.d | psxy -R -J -O -K -SE -Wthinnest >> ../example_18.ps

# Then draw 10 mGal contours and overlay 50 mGal contour in green

grdcontour AK_gulf_grav.nc -J -C20 -B2f1WSEn -O -K -U/-1.25i/-0.75i/"Example 18 in Cookbook" \
	-Y-4.85i >> ../example_18.ps
grdcontour AK_gulf_grav.nc -J -C10 -L49/51 -O -K -Dsm -Wcthin,green >> ../example_18.ps
pscoast -R -J -O -K -Di -Ggray -Wthinnest >> ../example_18.ps
$AWK '{print $1, $2, 0, 400, 400}' pratt.d | psxy -R -J -O -K -SE -Wthinnest >> ../example_18.ps
\rm -f sm_*[0-9].xyz	# Only consider closed contours

# Now determine centers of each enclosed seamount > 50 mGal but only plot
# the ones within 200 km of Pratt seamount.

# First determine mean location of each closed contour and
# add it to the file centers.d

\rm -f centers.d
foreach file (sm_*.xyz)
	$AWK 'BEGIN{x=0;y=0;n=0};{x+=$1;y+=$2;n++};END{print x/n,y/n}' $file >>! centers.d
end

# Only plot the ones within 200 km

gmtselect -C200/pratt.d centers.d -fg | \
   psxy -R -J -O -K -SC0.04i -Gred -Wthinnest >> ../example_18.ps
psxy -R -J -O -K -ST0.1i -Gyellow -Wthinnest pratt.d >> ../example_18.ps

# Then report the volume and area of these seamounts only
# by masking out data outside the 200 km-radius circle
# and then evaluate area/volume for the 50 mGal contour

grdmath -R `cat pratt.d` GDIST = mask.nc
grdclip mask.nc -Sa200/NaN -Sb200/1 -Gmask.nc
grdmath AK_gulf_grav.nc mask.nc MUL = tmp.nc
set info = `grdvolume tmp.nc -C50 -Sk`

psxy -R -J -A -O -K -L -Wthin -Gwhite << END >> ../example_18.ps
-148.5	52.75
-140.5	52.75
-140.5	53.75
-148.5	53.75
END
pstext -R -J -O << END >> ../example_18.ps
-148 53.08 14 0 1 LM Areas: $info[2] km@+2@+
-148 53.42 14 0 1 LM Volumes: $info[3] mGal\264km@+2@+
END

# Clean up

\rm -f grav.cpt sm_*.xyz *_i.nc tmp.nc mask.nc pratt.d center* .gmt*
