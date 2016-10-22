#!/bin/bash
#	$Id: GMT_tut_19.sh 16573 2016-06-19 02:42:45Z pwessel $
#
gmt makecpt -Ctopo -T1000/5000 > t.cpt
gmt grdgradient "${tut:-../tutorial}"/us.nc -Ne0.8 -A100 -fg -Gus_i.nc
gmt grdview "${tut:-../tutorial}"/us.nc -JM6i -p135/35 -Qi50 -Ius_i.nc -Ct.cpt -Ba -JZ0.5i > GMT_tut_19.ps
