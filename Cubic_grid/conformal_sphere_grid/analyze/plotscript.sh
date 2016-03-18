#/bin/bash

gnuplot <<EOF
set term png
set output "angle.png"
plot "angle_distribution.dat" using 1:2 with impulses ti "Angle distribution"
EOF

gnuplot <<EOF
set term png
set output "distance.png"
plot "distance_distribution.dat" using 1:2 with impulses ti "Linear size distribution"
EOF

gnuplot <<EOF
set term png
set output "square.png"
plot "square_distribution.dat" using 1:2 with impulses ti "Square size distribution"
EOF

gnuplot <<EOF
set xrange [89:91]
set term png
set output "angle_zoom.png"
plot "./angle_distribution.dat" using 1:2 with impulses ti "Angle distribution"
EOF

