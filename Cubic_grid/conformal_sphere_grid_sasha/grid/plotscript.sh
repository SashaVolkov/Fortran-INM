#/bin/bash

gnuplot <<EOF
set term png size 1280,960
set output "edges.png"
set hidden3d

set view 60,60,1.5,1
set xrange[-1:1]
set yrange[-1:1]
set zrange[-1:1]
set ticslevel 0
splot "edge1.dat" using 1:2:3 with lines noti,"edge2.dat" using 1:2:3 with lines noti,"edge3.dat" using 1:2:3 with lines noti,"edge4.dat" using 1:2:3 with lines noti,"edge5.dat" using 1:2:3 with lines noti,"edge6.dat" using 1:2:3 with lines noti

set output "zoom_edges.png"
set hidden3d
set view 60,60,1.5,1

set ticslevel 0
set xrange[0.5:0.7]
set yrange[0.5:0.7]
set zrange[0:1]

splot "edge1.dat" using 1:2:3 with lines noti,"edge2.dat" using 1:2:3 with lines noti,"edge3.dat" using 1:2:3 with lines noti,"edge4.dat" using 1:2:3 with lines noti,"edge5.dat" using 1:2:3 with lines noti,"edge6.dat" using 1:2:3 with lines noti

EOF
