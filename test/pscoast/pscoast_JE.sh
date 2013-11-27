#!/bin/bash
#
#	$Id: pscoast_JE.sh 12114 2013-09-03 19:19:00Z fwobbe $
# Make sure when fixed it works for all resolutions -D?

ps=pscoast_JE.ps

gmt pscoast -JE13:25/52:31/10/7i -Rg -Gred -Sblue -Dl -P > $ps

