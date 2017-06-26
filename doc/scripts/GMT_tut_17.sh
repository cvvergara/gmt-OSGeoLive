#!/bin/bash
#	$Id: GMT_tut_17.sh 18240 2017-05-27 23:39:17Z pwessel $
#
gmt makecpt -Cno_green -T-2/30/2 > otemp.cpt
gmt grdimage -Rg -JW180/9i "@otemp.anal1deg.nc?otemp[2,0]" -Cotemp.cpt -Bag > GMT_tut_17.ps
