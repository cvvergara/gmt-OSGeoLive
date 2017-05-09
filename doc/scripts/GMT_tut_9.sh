#!/bin/bash
#	$Id: GMT_tut_9.sh 17436 2017-01-13 00:22:07Z pwessel $
#
gmt pscoast -R130/150/35/50 -JM6i -B5 -P -Ggray -K > GMT_tut_9.ps
gmt psxy -R -J -O "${tut:-../tutorial}"/quakes.ngdc -Wfaint -i4,3,5,6+s0.1 -h3 -Scc -C"${tut:-../tutorial}"/quakes.cpt >> GMT_tut_9.ps
