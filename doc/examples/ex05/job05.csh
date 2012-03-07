#!/bin/csh
#		$Id: job05.csh,v 1.7 2011/03/01 01:34:48 remko Exp $
#		GMT EXAMPLE 05
#
# Purpose:	Generate grid and show monochrome 3-D perspective
# GMT progs:	grdgradient, grdmath, grdview, pstext
# Unix progs:	echo, rm

grdmath -R-15/15/-15/15 -I0.3 X Y HYPOT DUP 2 MUL PI MUL 8 DIV COS EXCH NEG 10 DIV EXP MUL \
	= sombrero.nc
echo '-5	128	5	128' >! gray.cpt
grdgradient sombrero.nc -A225 -Gintensity.nc -Nt0.75
grdview sombrero.nc -JX6i -JZ2i -B5/5/0.5SEwnZ -N-1/white -Qs -Iintensity.nc -X1.5i -K \
	-Cgray.cpt -R-15/15/-15/15/-1/1 -E120/30 -U/-1.25i/-0.75i/"Example 5 in Cookbook" \
	>! ../example_05.ps
echo "4.1 5.5 50 0 33 BC z(r) = cos (2@~p@~r/8) * e@+-r/10@+" | pstext -R0/11/0/8.5 -Jx1i -O \
	>> ../example_05.ps

\rm -f gray.cpt sombrero.nc intensity.nc .gmt*
