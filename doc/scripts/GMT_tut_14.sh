#!/bin/bash
#	$Id: GMT_tut_14.sh 16573 2016-06-19 02:42:45Z pwessel $
#
gmt makecpt -Crainbow -T-20/60/10 > disc.cpt
gmt makecpt -Crainbow -T-20/60 > cont.cpt
gmt psbasemap -R0/6/0/9 -Jx1i -P -B0 -K -Xc > GMT_tut_14.ps
gmt psscale -Dx1i/1i+w4i/0.5i+h -Cdisc.cpt -Ba -B+tdiscrete -O -K >> GMT_tut_14.ps
gmt psscale -Dx1i/3i+w4i/0.5i+h -Ccont.cpt -Ba -B+tcontinuous -O -K >> GMT_tut_14.ps
gmt psscale -Dx1i/5i+w4i/0.5i+h -Cdisc.cpt -Ba -B+tdiscrete -I0.5 -O -K >> GMT_tut_14.ps
gmt psscale -Dx1i/7i+w4i/0.5i+h -Ccont.cpt -Ba -B+tcontinuous -I0.5 -O >> GMT_tut_14.ps
