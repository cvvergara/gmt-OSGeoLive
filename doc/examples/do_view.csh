#!/bin/csh
#
#	$Id: do_view.csh 9545 2011-07-27 19:31:54Z pwessel $
#
#	Simple driver to view all examples using ghostview
#
if ($#argv == 1) then
	set viewer = $1
else
	set viewer = gv
endif

foreach f (example_*.ps)
	$viewer $f
end
