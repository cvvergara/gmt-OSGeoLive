#!/bin/csh
#		$Id: job09.csh 9545 2011-07-27 19:31:54Z pwessel $
#		GMT EXAMPLE 09
#
# Purpose:	Make wiggle plot along track from geoid deflections
# GMT progs:	pswiggle, pstext, psxy
# Unix progs:	$AWK, ls, paste, tail, rm

pswiggle track_*.xys -R185/250/-68/-42 -U"Example 9 in Cookbook" -K -Jm0.13i -Ba10f5 -Gblack -Z2000 \
	-Wthinnest -S240/-67/500/@~m@~rad >! ../example_09.ps
psxy -R -J -O -K ridge.xy -Wthicker >> ../example_09.ps
psxy -R -J -O -K -m fz.xy -Wthinner,- >> ../example_09.ps
\rm -f tmp
# Make label file
foreach file (track_*.xys)
	tail -1 $file >> tmp
end
ls -1 track_*.xys | $AWK -F. '{print $2}' >! tracks.lis
paste tmp tracks.lis | $AWK '{print $1, $2, 10, 50, 1, "RM", $4}' | pstext -R -J -D-0.05i/-0.05i -O \
	>> ../example_09.ps

\rm -f tmp tracks.lis .gmt*
