#!/bin/bash
#	$Id: GMT_transverse_merc.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt pscoast -R20/30/50/45r -Jt35/0.18i -Bag -Dl -A250 -Glightbrown -Wthinnest -P \
	-Sseashell > GMT_transverse_merc.ps
