#!/bin/bash
#	$Id: GMT_tut_16.sh 16573 2016-06-19 02:42:45Z pwessel $
#
gmt makecpt -Crainbow -T1000/5000 > topo.cpt
gmt grdgradient "${tut:-../tutorial}"/us.nc -Ne0.8 -A100 -fg -Gus_i.nc
gmt grdimage "${tut:-../tutorial}"/us.nc -Ius_i.nc -JM6i -P -Ba -Ctopo.cpt -K > GMT_tut_16.ps
gmt psscale -DjTC+w5i/0.25i+h+o0/-1i -R"${tut:-../tutorial}"/us.nc -J -Ctopo.cpt -I0.4 -Bxa -By+lm -O >> GMT_tut_16.ps
