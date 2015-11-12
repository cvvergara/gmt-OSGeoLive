#!/bin/bash
#	$Id: GMT_tut_2.sh 14523 2015-07-15 02:26:00Z pwessel $
#
gmt psbasemap -R1/10000/1e20/1e25 -JX9il/6il  -Bxa2+l"Wavelength (m)" -Bya1pf3+l"Power (W)" -BWS > GMT_tut_2.ps
