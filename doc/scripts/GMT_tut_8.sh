#!/bin/bash
#	$Id: GMT_tut_8.sh 14538 2015-07-15 19:43:07Z pwessel $
#
gmt psxy "${tut:-../tutorial}"/data -R0/6/0/6 -Jx1i -Baf -P -K -Wthinner > GMT_tut_8.ps
gmt psxy "${tut:-../tutorial}"/data -R -J -O -W -Si0.2i >> GMT_tut_8.ps