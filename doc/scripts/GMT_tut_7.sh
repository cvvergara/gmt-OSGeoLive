#!/bin/bash
#	$Id: GMT_tut_7.sh 14538 2015-07-15 19:43:07Z pwessel $
#
gmt psxy "${tut:-../tutorial}"/data -R0/6/0/6 -Jx1i -P -Baf -Wthinner > GMT_tut_7.ps
