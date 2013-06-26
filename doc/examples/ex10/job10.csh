#!/bin/csh
#		$Id: job10.csh 9545 2011-07-27 19:31:54Z pwessel $
#		GMT EXAMPLE 10
#
# Purpose:	Make 3-D bar graph on top of perspective map
# GMT progs:	pscoast, pstext, psxyz
# Unix progs:	$AWK, rm

pscoast -Rd -JX8id/5id -Dc -Gblack -E200/40 -K -U"Example 10 in Cookbook" \
	>! ../example_10.ps
psxyz agu2008.d -R-180/180/-90/90/1/100000 -J -JZ2.5il -So0.3ib1 -Ggray -Wthinner -O -K -E200/40 \
	-B60g60/30g30/a1p:Memberships:WSneZ >> ../example_10.ps
$AWK '{print $1, $2, 20, 0, 0, "RM", $3}' agu2008.d | pstext -Rd -J -O -K -E200/40 \
	-Gwhite -Sthinner -D-0.2i/0 >> ../example_10.ps
echo "4.5 6 30 0 5 BC AGU 2008 Membership Distribution" | pstext -R0/11/0/8.5 -Jx1i -O \
	>> ../example_10.ps

\rm -f .gmt*
