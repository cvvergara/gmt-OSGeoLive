#!/bin/sh
#
#	$Id: do_view.sh,v 1.5 2011/06/24 01:23:49 guru Exp $
#
#	Simple driver to view all examples using ghostview
#
viewer=${1:-gv}
for f in example_*.ps
do
	$viewer $f
done
