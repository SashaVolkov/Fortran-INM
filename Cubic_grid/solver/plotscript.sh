#/bin/bash

line=$(cat init | sed 's/ //g')
args=(${line//,/ })
dim="${args[0]}"
grid_type="${args[5]}"
step="${args[6]}"
test_number="${args[7]}"

if [[ $grid_type == 1 ]]; then
	grid="equiang"
elif [ $grid_type == 0 ]&&[ $rescale == 0 ]; then
	grid="simple"
elif [ $grid_type == 0 ] && [ $rescale == 1 ]; then
	grid="tan"
fi

cd datFiles/"$(( 2*$dim ))"/"$(( 2*$step ))"th/"$grid"/test"$test_number"
pwd

gnuplot <<EOF
set term png size 800, 800
set output "../../../../pic/$(( 2*$step ))th/CFL_$(( 2*$dim )).png"
set xlabel "Days"
set key inside left top vertical Right noreverse
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "CFL.dat" w l ti "CFL_{$grid_type}"
EOF

gnuplot <<EOF
set term png size 800, 800
set output "../../../../pic/$(( 2*$step ))th/L_equiang_$(( 2*$dim ))_h.png"
set xlabel "Days"
set key inside left top vertical Right noreverse
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "L1_h.dat" w l ti "L_1", "L2_h.dat" w l ti "L_2", "C_h.dat" w l ti "C"
EOF


gnuplot <<EOF
set term png size 800, 800
set output "../../../../pic/$(( 2*$step ))th/L_equiang_$(( 2*$dim ))_v.png"
set xlabel "Days"
set key inside left top vertical Right noreverse
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "L1_v.dat" w l ti "L_1", "L2_v.dat" w l ti "L_2", "C_v.dat" w l ti "C"
EOF


gnuplot <<EOF
set term png size 800, 800
set output "../../../../pic/$(( 2*$step ))th/L1_$(( 2*$dim ))_v.png"
set xlabel "Days"
set key inside left top vertical Right noreverse
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "L1_v.dat" w l ti "L_1equiang"
EOF


gnuplot <<EOF
set term png size 800, 800
set output "../../../../pic/$(( 2*$step ))th/L1_$(( 2*$dim ))_h.png"
set xlabel "Days"
set key inside left top vertical Right noreverse
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "L1_h.dat" w l ti "L_1equiang"
EOF


gnuplot <<EOF
set term png size 800, 800
set output "../../../../pic/$(( 2*$step ))th/L2_$(( 2*$dim ))_h.png"
set key inside right top vertical Right noreverse
set xlabel "Days"
set ylabel "Normalized error"
set mxtics 5
set mytics 5
set grid ytics lt 2 lw 2 lc rgb "#bbbbbb"
set grid xtics lt 2 lw 2 lc rgb "#bbbbbb"
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "L2_h.dat" w l lw 1 ti " L_2 equiang"
EOF


gnuplot <<EOF
set term png size 800, 800
set output "../../../../pic/$(( 2*$step ))th/L2_$(( 2*$dim ))_v.png"
set key inside right top vertical Right noreverse
set xlabel "Days"
set ylabel "Normalized error"
set mxtics 5
set mytics 5
set grid ytics lt 2 lw 2 lc rgb "#bbbbbb"
set grid xtics lt 2 lw 2 lc rgb "#bbbbbb"
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "L2_v.dat" w l lw 1 ti " L_2 equiang"
EOF


gnuplot <<EOF
set term png size 800, 800
set output "../../../../pic/$(( 2*$step ))th/L_inf_$(( 2*$dim ))_h.png"
set key inside left top vertical Right noreverse
set xlabel "Days"
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
plot "C_h.dat" w l ti "C_{equiang}"
EOF



gnuplot <<EOF
set term png size 800, 800
set output "../../../../pic/$(( 2*$step ))th/L_inf_$(( 2*$dim ))_v.png"
set key inside left top vertical Right noreverse
set xlabel "Days"
set grid ytics lt 0 lw 1 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 1 lc rgb "#bbbbbb"
set mxtics 5
set mytics 5
plot "C_v.dat" w l ti "C_{equiang}"
EOF


# cd ../..

# gnuplot <<EOF
# set term png size 800, 800
# set output "pic/$(( 2*$step ))th/L2_equiang_h.png"
# set key inside left top vertical Right noreverse
# set xlabel "Days"
# set mxtics 5
# set mytics 5
# set grid ytics lt 2 lw 2 lc rgb "#bbbbbb"
# set grid xtics lt 2 lw 2 lc rgb "#bbbbbb"
# set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
# set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
# plot "40/$(( 2*$step ))th/equiang/L2_h.dat" w l ti " L_2 40", "60/$(( 2*$step ))th/equiang/L2_h.dat" w l ti " L_2 60", "80/$(( 2*$step ))th/equiang/L2_h.dat" w l ti " L_2 80"
# EOF


gnuplot <<EOF
set term png size 800, 800
set output "pic/$(( 2*$step ))th/L2_equiang_v.png"
set key inside left top vertical Right noreverse
set xlabel "Days"
set mxtics 5
set mytics 5
set grid ytics lt 2 lw 2 lc rgb "#bbbbbb"
set grid xtics lt 2 lw 2 lc rgb "#bbbbbb"
set grid mxtics lt 0 lw 1 lc rgb "#bbbbbb"
set grid mytics lt 0 lw 1 lc rgb "#bbbbbb"
plot "40/$(( 2*$step ))th/equiang/L2_v.dat" w l ti " L_2 40", "60/$(( 2*$step ))th/equiang/L2_v.dat" w l ti " L_2 60", "80/$(( 2*$step ))th/equiang/L2_v.dat" w l ti " L_2 80"
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