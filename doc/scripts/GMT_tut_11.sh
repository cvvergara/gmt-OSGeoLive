#!/bin/bash
#	$Id: GMT_tut_11.sh 14538 2015-07-15 19:43:07Z pwessel $
#
gmt grdcontour "${tut:-../tutorial}"/bermuda.nc -JM6i -C250 -A1000 -P -Ba > GMT_tut_11.ps
