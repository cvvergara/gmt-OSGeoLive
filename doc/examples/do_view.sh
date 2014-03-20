#!/bin/sh
#
#	$Id: do_view.sh 9545 2011-07-27 19:31:54Z pwessel $
#
#	Simple driver to view all examples using ghostview
#
viewer=${1:-gv}
for f in example_*.ps
do
	$viewer $f
done
