#/bin/bash

cd datFiles

gnuplot <<EOF
set term png
set output "CFL.png"
set xlabel "Days"
plot "CFL_conf_tan.dat" w l ti "CFL_tan", "CFL_conf_simple.dat" w l ti "CFL_simple", "CFL_equiang.dat" w l ti "CFL_equiang"
EOF

gnuplot <<EOF
set term png
set output "L_conf_simple.png"
set xlabel "Days"
plot "L1_conf_simple.dat" w l ti "L1", "L2_conf_simple.dat" w l ti "L2", "L_inf_conf_simple.dat" w l ti "L_inf"
EOF


gnuplot <<EOF
set term png
set output "L_conf_tan.png"
set xlabel "Days"
plot "L1_conf_tan.dat" w l ti "L1", "L2_conf_tan.dat" w l ti "L2", "L_inf_conf_tan.dat" w l ti "L_inf"
EOF

gnuplot <<EOF
set term png
set output "L_equiang.png"
set xlabel "Days"
plot "L1_equiang.dat" w l ti "L1", "L2_equiang.dat" w l ti "L2", "L_inf_equiang.dat" w l ti "L_inf"
EOF


gnuplot <<EOF
set term png
set output "L1.png"
set xlabel "Days"
plot "L1_conf_simple.dat" w l ti "L1_simple", "L1_conf_tan.dat" w l ti "L1_tan", "L1_equiang.dat" w l ti "L1_equiang"
EOF

gnuplot <<EOF
set term png
set output "L2.png"
set xlabel "Days"
plot "L2_conf_simple.dat" w l ti "L2_simple", "L2_conf_tan.dat" w l ti "L2_tan", "L2_equiang.dat" w l ti "L2_equiang"
EOF

gnuplot <<EOF
set term png
set output "L_inf.png"
set xlabel "Days"
plot "L_inf_conf_simple.dat" w l ti "L_inf_simple", "L_inf_conf_tan.dat" w l ti "L_inf_tan", "L_inf_equiang.dat" w l ti "L_inf_equiang"
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