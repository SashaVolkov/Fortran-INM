#/bin/bash

rm -rf *.o *.mod *.out *.file analyze *~ 2>/dev/null

Files="sphere_geometry.f90 diagn.f90 scan-print.f90 grid_interp.f90"
Files=$Files" Regridder.f90"


netcdf="/data4t/avolkov/util/netcdf-2016Jan-13.1"
netcdf="/home/sasha/netcdf"

 # -check all -traceback -ftrapuv
mpiifort -openmp -O3 $Files -I $netcdf/inc -L $netcdf/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm 2> err.file
# /home/sasha/Fortran/Comands/./compo geometry.o conformal.o matmul.o morphism.o grid_generator.o data_analyzer.o spherical.o main.o
	echo "compilation status" $?

if [[ `grep -c error err.file` > 0 ]]; then
	echo "Look for" `grep -c error err.file` "errors in err.file"
	echo `grep -c warning err.file` "warnings"
else

	echo `grep -c error err.file` "errors"
	echo `grep -c warning err.file` "warnings"

	if [[ $1 != "compile" ]]; then
		export OMP_NUM_THREADS=4
		time mpiexec -n 1 ./a.out

	fi


fi
