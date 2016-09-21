#/bin/bash

rm -rf *.o *.mod *.out *.file analyze *~ 2>/dev/null

Files="geometry.f90 matmul.f90 simple_rotations.f90 spherical.f90"
Files=$Files" projections.f90 matrix_rotation.f90"
Files=$Files" parallel_cubic.f90 grid_generator.f90 metrics.f90 grid_var.f90 derivatives.f90 interpolation.f90 func_var.f90 messenger.f90 printer.f90 schemes.f90 diagnostic.f90"
Files=$Files" Solver_shallow_water.f90"

netcdf="/data4t/avolkov/util/netcdf-2016Jan-13.1"
netcdf="/home/sasha/netcdf"

 # -check all -traceback -ftrapuv
mpiifort -O3 $Files -I $netcdf/inc -L $netcdf/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm 2> err.file
# /home/sasha/Fortran/Comands/./compo geometry.o conformal.o matmul.o morphism.o grid_generator.o data_analyzer.o spherical.o main.o
	echo "compilation status" $?

if [[ `grep -c error err.file` > 0 ]]; then
	echo "Look for" `grep -c error err.file` "errors in err.file"
	echo `grep -c warning err.file` "warnings"
else

	echo `grep -c error err.file` "errors"
	echo `grep -c warning err.file` "warnings"

	if [[ $1 != "compile" ]]; then
		export OMP_NUM_THREADS=1
		mpiexec -n $1 ./a.out

		./plotscript.sh
	fi

echo "Regridding"
cd regrid
./run_intel_reg.sh

fi
