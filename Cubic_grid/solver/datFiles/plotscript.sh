#/bin/bash

gnuplot <<EOF
set term png
set output "CFL_x.png"
set xlabel "Time"
plot "CFL_x_tan.dat" w l ti "CFL_x_tan", "CFL_x_simple.dat" w l ti "CFL_x_simple"
EOF

gnuplot <<EOF
set term png
set output "CFL_y.png"
set xlabel "Time"
plot "CFL_y_tan.dat" w l ti "CFL_y_tan", "CFL_y_simple.dat" w l ti "CFL_y_simple"
EOF



gnuplot <<EOF
set term png
set output "L1.png"
set xlabel "Time"
plot "L1_tan.dat" w l ti "L1_tan", "L1_simple.dat" w l ti "L1_simple"
EOF

gnuplot <<EOF
set term png
set output "L2.png"
set xlabel "Time"
plot "L2_tan.dat" w l ti "L2_tan", "L2_simple.dat" w l ti "L2_simple"
EOF

gnuplot <<EOF
set term png
set output "L_inf.png"
set xlabel "Time"
plot "L_inf_tan.dat" w l ti "L_inf_tan", "L_inf_simple.dat" w l ti "L_inf_simple"
EOF

gnuplot <<EOF
set term png
set output "angle_tan.png"
plot "angle_distribution_tan.dat" using 1:2 with impulses ti "Angle distribution"
EOF


gnuplot <<EOF
set term png
set output "angle_simple.png"
plot "angle_distribution_simple.dat" using 1:2 with impulses ti "Angle distribution"
EOF

gnuplot <<EOF
set xrange [89:91]
set term png
set output "angle_tan_zoom.png"
plot "./angle_distribution_tan.dat" using 1:2 with impulses ti "Angle distribution"
EOF

gnuplot <<EOF
set xrange [89:91]
set term png
set output "angle_simple_zoom.png"
plot "./angle_distribution_simple.dat" using 1:2 with impulses ti "Angle distribution"
EOF