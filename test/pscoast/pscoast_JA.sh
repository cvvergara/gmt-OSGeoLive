#!/bin/bash
#
#	$Id: pscoast_JA.sh 12114 2013-09-03 19:19:00Z fwobbe $
# Make sure when fixed it works for all resolutions -D?

ps=pscoast_JA.ps

gmt pscoast -JA13:25/52:31/10/7i -Rg -Gred -Sblue -Dl -P > $ps

