#!/bin/bash
#	$Id: GMT_obl_merc.sh 11490 2013-05-16 06:26:21Z pwessel $
#
gmt pscoast -R270/20/305/25r -JOc280/25.5/22/69/4.8i -Bag -Di -A250 -Gburlywood -Wthinnest -P \
	-Tf301.5/23/0.4i/2 -Sazure --FONT_TITLE=8p --MAP_TITLE_OFFSET=0.05i > GMT_obl_merc.ps