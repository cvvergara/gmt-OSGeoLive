#!/bin/bash
#	$Id: GMT_tut_11.sh 18314 2017-06-03 03:26:22Z pwessel $
#
gmt grdcontour @tut_bathy.nc -JM6i -C250 -A1000 -P -Ba > GMT_tut_11.ps
