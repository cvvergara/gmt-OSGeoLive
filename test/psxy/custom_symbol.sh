#!/bin/bash
#       $Id: custom_symbol.sh 12350 2013-10-17 13:38:22Z fwobbe $
#
# Check two custom symbols symbols with new variables and text capabilities

ps=custom_symbol.ps

sed -n 1p "${src:-.}"/custom_data.txt | gmt psxy -R0/20/10/30 -JM6i -P -Bag -W0.5p,red -Sk"${src:-.}"/comet/2.5i -K -Xc -Yc > $ps
sed -n 2p "${src:-.}"/custom_data.txt | gmt psxy -R -J -W1p,green -Sk"${src:-.}"/dip/2.5i -O >> $ps
