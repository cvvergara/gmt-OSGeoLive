#!/bin/bash
#	$Id: GMT_tut_1.sh 14523 2015-07-15 02:26:00Z pwessel $
#
gmt psbasemap -R10/70/-3/8 -JX4i/3i -Ba -B+glightred+t"My first plot" -P > GMT_tut_1.ps
