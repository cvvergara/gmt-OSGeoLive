#!/bin/bash
#
#	$Id: withgdal.sh 12350 2013-10-17 13:38:22Z fwobbe $

GDAL=`gmt grdreformat 2>&1 | grep -c gd`
if [ $GDAL -eq 0 ]; then exit; fi
	
ps=withgdal.ps

# RGB image
gmt psimage "${src:-.}"/../grdimage/gdal/needle.jpg -W7c -P -Y15c -K > $ps

# Same image as above but as idexed
gmt psimage "${src:-.}"/../grdimage/gdal/needle.png -W7c -X7.5c -O -K >> $ps

# Convert RGB to YIQ
gmt psimage "${src:-.}"/../grdimage/gdal/needle.jpg -M -W7c -X-7.5c -Y-7c -O -K >> $ps

# Convert Indexed to YIQ
gmt psimage "${src:-.}"/../grdimage/gdal/needle.png -M -W7c -X7.5c -O -K >> $ps

# A gray image (one band, no color map)
gmt psimage "${src:-.}"/../grdimage/gdal/vader.jpg -W4c -X-2.5c -Y4.5c -O >> $ps

