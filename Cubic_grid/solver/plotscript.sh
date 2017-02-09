#/bin/bash

line=$(cat init | sed 's/ //g')
args=(${line//,/ })
dim="${args[0]}"
step="${args[6]}"

cd datFiles/"$(( 2*$dim ))"/"$(( 2*$step ))"th
pwd

gnuplot <<EOF
set term png size 800, 800
set output "../../pic/$(( 2*$step ))th/CFL_$(( 2*$dim )).png"
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
set output "../../pic/$(( 2*$step ))th/L_conf_simple_$(( 2*$dim )).png"
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
set output "../../pic/$(( 2*$step ))th/L_conf_tan_$(( 2*$dim )).png"
set xlabel "Cycles"
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
set output "../../pic/$(( 2*$step ))th/L_equiang_$(( 2*$dim )).png"
set xlabel "Cycles"
set key inside left top vertical Right noreverse
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set yrange [0:.2]
set mxtics 5
set mytics 5
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "equiang/L1.dat" w l ti "L1", "equiang/L2.dat" w l ti "L2", "equiang/C.dat" w l ti "C"
EOF


gnuplot <<EOF
set term png size 800, 800
set output "../../pic/$(( 2*$step ))th/L1_$(( 2*$dim )).png"
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

# set terminal postscript eps size 4, 3
gnuplot <<EOF
set term png size 800, 800
set output "../../pic/$(( 2*$step ))th/L2_$(( 2*$dim )).png"
set key inside right top vertical Right noreverse
set xlabel "Cycles"
set ylabel "Normalized error"
set yrange [0:.1]
set xrange [0:10]
set mxtics 5
set mytics 5
set grid ytics lt 2 lw 2 lc rgb "#bbbbbb"
set grid xtics lt 2 lw 2 lc rgb "#bbbbbb"
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "simple/L2.dat" w l lw 5 ti "L2 conf", "tan/L2.dat" w l lw 5 lt 5 ti "L2 tan", "equiang/L2.dat" w l lw 5 lt 3 ti "L2 equiang"
EOF

gnuplot <<EOF
set term png size 800, 800
set output "../../pic/$(( 2*$step ))th/L_inf_$(( 2*$dim )).png"
set key inside left top vertical Right noreverse
set xlabel "Cycles"
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
plot "simple/C.dat" w l ti "C_simple", "tan/C.dat" w l ti "C_tan", "equiang/C.dat" w l ti "C_equiang"
EOF


cd ../..

gnuplot <<EOF
set term png size 800, 800
set output "pic/$(( 2*$step ))th/L2_equiang.png"
set key inside left top vertical Right noreverse
set xlabel "Cycles"
set mxtics 5
set mytics 5
set grid ytics lt 2 lw 2 lc rgb "#bbbbbb"
set grid xtics lt 2 lw 2 lc rgb "#bbbbbb"
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "40/$(( 2*$step ))th/equiang/L2.dat" w l ti "L2_40", "60/$(( 2*$step ))th/equiang/L2.dat" w l ti "L2_60", "80/$(( 2*$step ))th/equiang/L2.dat" w l ti "L2_80"
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