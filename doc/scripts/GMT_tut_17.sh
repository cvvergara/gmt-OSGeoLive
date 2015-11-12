#!/bin/bash
#	$Id: GMT_tut_17.sh 14538 2015-07-15 19:43:07Z pwessel $
#
gmt makecpt -Cno_green -T-2/30/2 > otemp.cpt
gmt grdimage -Rg -JW180/9i "${tut:-../tutorial}/otemp.anal1deg.nc?otemp[2,0]" -Cotemp.cpt -Bag > GMT_tut_17.ps
