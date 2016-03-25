#/bin/bash

gnuplot <<EOF
set term png size 1280,960
set output "faces.png"
set hidden3d

set view 60,60,1.5,1
set xrange[-1:1]
set yrange[-1:1]
set zrange[-1:1]
set ticslevel 0
splot "face1.dat" using 1:2:3 with lines noti,"face2.dat" using 1:2:3 with lines noti,"face3.dat" using 1:2:3 with lines noti,"face4.dat" using 1:2:3 with lines noti,"face5.dat" using 1:2:3 with lines noti,"face6.dat" using 1:2:3 with lines noti

set output "zoom_faces.png"
set hidden3d
set view 60,60,1.5,1

set ticslevel 0
set xrange[0.5:0.7]
set yrange[0.5:0.7]
set zrange[0:1]

splot "face1.dat" using 1:2:3 with lines noti,"face2.dat" using 1:2:3 with lines noti,"face3.dat" using 1:2:3 with lines noti,"face4.dat" using 1:2:3 with lines noti,"face5.dat" using 1:2:3 with lines noti,"face6.dat" using 1:2:3 with lines noti

EOF
