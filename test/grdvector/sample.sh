#!/bin/bash
# $Id: sample.sh 12480 2013-11-09 17:02:40Z fwobbe $
# Test sub-sampling of grids for vectors
ps=sample.ps
gmt pscoast -R60/105/-20/40 -JM95.0/35/16c -Gbisque -K -Bafg8 -P -Xc > $ps
gmt grdvector nuvel1.vx.1.5.nc nuvel1.vy.1.5.nc -T -R -J -I2 -Si150.0 -Q0.2i+e -Wthicker,lightgray -Glightgray -O -K -t50 >> $ps
gmt grdvector nuvel1.vx.1.5.nc nuvel1.vy.1.5.nc -T -R -J -I8 -Si150.0 -Q0.2i+e -Wthicker,black -Gred -O >> $ps
