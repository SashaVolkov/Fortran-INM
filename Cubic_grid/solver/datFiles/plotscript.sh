#/bin/bash

gnuplot <<EOF
set term png
set output "CFL.png"
set xlabel "Days"
plot "CFL_tan.dat" w l ti "CFL_tan", "CFL_simple.dat" w l ti "CFL_simple"
EOF

gnuplot <<EOF
set term png
set output "L.png"
set xlabel "Days"
plot "L1_simple.dat" w l ti "L1", "L2_simple.dat" w l ti "L2", "L_inf_simple.dat" w l ti "L_inf"
EOF

gnuplot <<EOF
set term png
set output "L1.png"
set xlabel "Days"
plot "L1_simple.dat" w l ti "L1_simple"
EOF

gnuplot <<EOF
set term png
set output "L2.png"
set xlabel "Days"
plot "L2_simple.dat" w l ti "L2_simple"
EOF

gnuplot <<EOF
set term png
set output "L_inf.png"
set xlabel "Days"
plot "L_inf_simple.dat" w l ti "L_inf_simple"
EOF

# gnuplot <<EOF
# set term png
# set output "angle_tan.png"
# set xlabel "Cell's angles Degrees"
# set ylabel "N"
# plot "angle_distribution_tan.dat" using 1:2 with impulses ti "Angle distribution"
# EOF


# gnuplot <<EOF
# set term png
# set output "angle_simple.png"
# set xlabel "Cell's angles Degrees"
# set ylabel "N"
# plot "angle_distribution_simple.dat" using 1:2 with impulses ti "Angle distribution"
# EOF

# gnuplot <<EOF
# set xrange [89:91]
# set term png
# set output "angle_tan_zoom.png"
# set xlabel "Cell's angles Degrees"
# set ylabel "N"
# plot "angle_distribution_tan.dat" using 1:2 with impulses ti "Angle distribution"
# EOF

# gnuplot <<EOF
# set xrange [89:91]
# set term png
# set output "angle_simple_zoom.png"
# set xlabel "Cell's angles Degrees"
# set ylabel "N"
# plot "angle_distribution_simple.dat" using 1:2 with impulses ti "Angle distribution"
# EOF

# gnuplot <<EOF
# set term png
# set output "cell_simple.png"
# set ylabel "N"
# set xlabel "Cell's area m^2"
# plot "cell_distribution_simple.dat" using 1:2 with impulses ti "Area distribution"
# EOF

# gnuplot <<EOF
# set term png
# set output "cell_tan.png"
# set ylabel "N"
# set xlabel "Cell's area m^2"
# plot "cell_distribution_tan.dat" using 1:2 with impulses ti "Area distribution"
# EOF


# gnuplot <<EOF
# set term png
# set output "dist_simple.png"
# set ylabel "N"
# set xlabel "Cell's side length m"
# plot "dist_distribution_simple.dat" using 1:2 with impulses ti "Distance distribution"
# EOF

# gnuplot <<EOF
# set term png
# set output "dist_tan.png"
# set ylabel "N"
# set xlabel "Cell's side length m"
# plot "dist_distribution_tan.dat" using 1:2 with impulses ti "Distance distribution"
# EOF