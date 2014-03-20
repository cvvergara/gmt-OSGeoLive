#!/bin/bash
#	$Id: GMT_pow.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt psxy -R0/100/0/10 -Jx0.3ip0.5/0.15i -Bxa1p -Bya2f1 -BWSne+givory -Wthick -P -K sqrt.d > GMT_pow.ps
gmt psxy -R -J -Sc0.075i -Ggreen -W -O sqrt.d10 >> GMT_pow.ps
