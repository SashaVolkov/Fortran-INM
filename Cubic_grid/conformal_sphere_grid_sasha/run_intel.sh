#/bin/bash

./clean.sh

mkdir grid 2>/dev/null
mkdir analyze 2>/dev/null

echo '#/bin/bash

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
' > analyze/plotscript.sh
chmod 755 analyze/plotscript.sh

echo '#/bin/bash

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

EOF' > grid/plotscript.sh
chmod 755 grid/plotscript.sh


mpiifort geometry.f90 mapping.f90 matmul.f90 simple_rotations.f90 projections.f90 matrix_rotation.f90 grid_generator.f90 data_analyzer.f90 spherical.f90 main.f90 2> err.file
# /home/sasha/Fortran/Comands/./compo geometry.o conformal.o matmul.o morphism.o grid_generator.o data_analyzer.o spherical.o main.o


if [[ `grep -c error err.file` > 0 ]]; then
	echo "Look for" `grep -c error err.file` "errors in err.file"
	echo `grep -c warning err.file` "warnings"
else

	echo `grep -c error err.file` "errors"
	echo `grep -c warning err.file` "warnings"

	mpiexec -n 1 ./a.out

	cd analyze
	./plotscript.sh
	cd ../grid
	./plotscript.sh
fi


