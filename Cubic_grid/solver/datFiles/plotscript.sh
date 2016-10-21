#/bin/bash

gnuplot <<EOF
set term png
set output "pic/CFL.png"
set xlabel "Days"
set xrange [0:6]
plot "simple/CFL.dat" w l ti "CFL_simple", "tan/CFL.dat" w l ti "CFL_tan", "equiang/CFL.dat" w l ti "CFL_equiang"
EOF

# gnuplot <<EOF
# set term png
# set output "CFL.png"
# set xlabel "Days"
# "CFL_conf_simple.dat" w l ti "CFL_simple", "CFL_equiang.dat" w l ti "CFL_equiang"
# EOF

# gnuplot <<EOF
# set term png
# set output "CFL.png"
# set xlabel "Days"
# "CFL_equiang.dat" w l ti "CFL_equiang"
# EOF

gnuplot <<EOF
set term png
set output "pic/L_conf_simple.png"
set xlabel "Days"
plot "simple/L1.dat" w l ti "L1", "simple/L2.dat" w l ti "L2", "simple/L_inf.dat" w l ti "L_inf"
EOF


gnuplot <<EOF
set term png
set output "pic/L_conf_tan.png"
set xlabel "Days"
plot "tan/L1.dat" w l ti "L1", "tan/L2.dat" w l ti "L2", "tan/L_inf.dat" w l ti "L_inf"
EOF

gnuplot <<EOF
set term png
set output "pic/L_equiang.png"
set xlabel "Days"
set xrange [0:10]
set yrange [0:18]
plot "equiang/L1.dat" w l ti "L1", "equiang/L2.dat" w l ti "L2", "equiang/L_inf.dat" w l ti "L_inf"
EOF


gnuplot <<EOF
set term png
set output "pic/L1.png"
set xlabel "Days"
plot "simple/L1.dat" w l ti "L1_simple", "tan/L1.dat" w l ti "L1_tan", "equiang/L1.dat" w l ti "L1_equiang"
EOF

gnuplot <<EOF
set term png
set output "pic/L2.png"
set xlabel "Days"
plot "simple/L2.dat" w l ti "L2_simple", "tan/L2.dat" w l ti "L2_tan", "equiang/L2.dat" w l ti "L2_equiang"
EOF

gnuplot <<EOF
set term png
set output "pic/L2eq.png"
set xlabel "Days"
set xrange [0:6]
plot "equiang/L2_2_2.dat" w l ti "L2_2_only", "equiang/L2.dat" w l ti "L2_2and4", "equiang/L2_4.dat" w l ti "L2_4_only"
EOF

gnuplot <<EOF
set term png
set output "pic/L_inf.png"
set xlabel "Days"
plot "simple/L_inf.dat" w l ti "L_inf_simple", "tan/L_inf.dat" w l ti "L_inf_tan", "equiang/L_inf.dat" w l ti "L_inf_equiang"
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