.. _psxy_ex_-L:

  ::

   #!/bin/bash
   ps=filler.ps
   cat << EOF > t.txt
   1 1
   2 3
   3 2
   4 4
   EOF
   psxy -R0/5/0/5 -JX3i -P -K -B0 t.txt -Gred -W2p -L+yb > $ps
   psxy -R -J -O -K -B0 t.txt -Gred -W2p -L+yt -X3.25i >> $ps
   psxy -R -J -O -K -B0 t.txt -Gred -W2p -L+xl -X-3.25i -Y3.25i >> $ps
   psxy -R -J -O -K -B0 t.txt -Gred -W2p -L+xr -X3.25i >> $ps
   psxy -R -J -O -K -B0 t.txt -Gred -W2p -L+y4 -X-3.25i -Y3.25i >> $ps
   psxy -R -J -O -K -B0 t.txt -Gred -W2p -L+x4.5 -X3.25i >> $ps
   psxy -R -J -O -T >> $ps
