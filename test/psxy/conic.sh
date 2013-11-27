#!/bin/bash
#	$Id: conic.sh 12114 2013-09-03 19:19:00Z fwobbe $
#
# Check clipping of line for a global conic plot

ps=conic.ps

gmt pscoast -R0/360/30/70 -JL180/50/40/60/6i -Gred -Dc -B30g30 -P -K > $ps
gmt psxy -R -J -O -W2p << EOF >> $ps
30	50
-30	45
EOF
