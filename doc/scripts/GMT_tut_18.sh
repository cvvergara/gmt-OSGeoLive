#!/bin/bash
#	$Id: GMT_tut_18.sh 18314 2017-06-03 03:26:22Z pwessel $
#
gmt grd2cpt @tut_bathy.nc -Cocean > bermuda.cpt
gmt grdview tut_bathy.nc -JM5i -P -JZ2i -p135/30 -Ba -Cbermuda.cpt > GMT_tut_18.ps
