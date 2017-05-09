#!/bin/bash
#	$Id: GMT_tut_15.sh 17806 2017-03-27 22:18:01Z pwessel $
#
gmt makecpt -Crainbow -T1000/5000 > topo.cpt
gmt grdimage "${tut:-../tutorial}"/us.nc -JM6i -P -Ba -BWSnE -Ctopo.cpt -K > GMT_tut_15.ps
gmt psscale -DJTC -R"${tut:-../tutorial}"/us.nc -J -Ctopo.cpt -I0.4 -Bxa -By+lm -O >> GMT_tut_15.ps
