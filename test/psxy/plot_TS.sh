#!/bin/bash
#	$Id: plot_TS.sh 12114 2013-09-03 19:19:00Z fwobbe $
#
# Plot lines with variable number of NaNs

ps=plot_TS.ps

gmt psxy chkPts_tseries.dat -i0,12 -JX24c/16c -Bxa24f3 -Bya2f1 -BWS -R0/324/12/26 -W1 -K > $ps
gmt psxy chkPts_tseries.dat -i0,11 -JX -R -W1,red -O >> $ps

