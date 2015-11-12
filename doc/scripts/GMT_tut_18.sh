#!/bin/bash
#	$Id: GMT_tut_18.sh 14538 2015-07-15 19:43:07Z pwessel $
#
gmt grd2cpt "${tut:-../tutorial}"/bermuda.nc -Cocean > bermuda.cpt
gmt grdview "${tut:-../tutorial}"/bermuda.nc -JM5i -P -JZ2i -p135/30 -Ba -Cbermuda.cpt > GMT_tut_18.ps
