#!/bin/bash
#	$Id: GMT_tut_12.sh 15178 2015-11-06 10:45:03Z fwobbe $
#
gmt nearneighbor -R245/255/20/30 -I5m -S40k -Gship.nc "${tut:-../tutorial}"/ship.xyz
gmt grdcontour ship.nc -JM6i -P -Ba -C250 -A1000 > GMT_tut_12.ps
