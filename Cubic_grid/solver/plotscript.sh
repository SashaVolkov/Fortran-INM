#/bin/bash

line=$(cat init | sed 's/ //g')
args=(${line//,/ })
dim="${args[0]}"
step="${args[6]}"

cd datFiles/"$(( 2*$dim ))"/"$(( 2*$step ))"th
pwd

gnuplot <<EOF
set term png size 800, 800
set output "pic/CFL.png"
set xlabel "Cycles"
set key inside left top vertical Right noreverse
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "tan/CFL.dat" w l ti "CFL_tan", "simple/CFL.dat" w l ti "CFL_simple", "equiang/CFL.dat" w l ti "CFL_equiang"
EOF

gnuplot <<EOF
set term png size 800, 800
set output "pic/L_conf_simple.png"
set yrange [0:2]
set xlabel "Cycles"
set key inside left top vertical Right noreverse
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "simple/L1.dat" w l ti "L1", "simple/L2.dat" w l ti "L2", "simple/C.dat" w l ti "C"
EOF


gnuplot <<EOF
set term png size 800, 800
set output "pic/L_conf_tan.png"
set xlabel "Cycles"
set yrange [0:2]
set key inside left top vertical Right noreverse
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "tan/L1.dat" w l ti "L1", "tan/L2.dat" w l ti "L2", "tan/C.dat" w l ti "C"
EOF

gnuplot <<EOF
set term png size 800, 800
set output "pic/L_equiang.png"
set xlabel "Cycles"
set key inside left top vertical Right noreverse
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "equiang/L1.dat" w l ti "L1", "equiang/L2.dat" w l ti "L2", "equiang/C.dat" w l ti "C"
EOF


gnuplot <<EOF
set term png size 800, 800
set output "pic/L1.png"
set xlabel "Cycles"
set key inside left top vertical Right noreverse
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "simple/L1.dat" w l ti "L1_simple", "tan/L1.dat" w l ti "L1_tan", "equiang/L1.dat" w l ti "L1_equiang"
EOF

gnuplot <<EOF
set term png size 800, 800
set output "pic/L2.png"
set key inside left top vertical Right noreverse
set xlabel "Cycles"
set mxtics 5
set mytics 5
set grid ytics lt 2 lw 2 lc rgb "#bbbbbb"
set grid xtics lt 2 lw 2 lc rgb "#bbbbbb"
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "simple/L2.dat" w l ti "L2_simple", "tan/L2.dat" w l ti "L2_tan", "equiang/L2.dat" w l ti "L2_equiang"
EOF

gnuplot <<EOF
set term png size 800, 800
set output "pic/L_inf.png"
set key inside left top vertical Right noreverse
set xlabel "Cycles"
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
plot "simple/C.dat" w l ti "C_simple", "tan/C.dat" w l ti "C_tan", "equiang/C.dat" w l ti "C_equiang"
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