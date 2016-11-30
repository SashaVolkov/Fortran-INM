#/bin/bash

line=$(cat init | sed 's/ //g')
args=(${line//,/ })
dim="${args[0]}"
rescale="${args[4]}"
grid_type="${args[5]}"

if [ ! -d mod_files ]; then
	mkdir mod_files
fi

rm -rf *.o mod_files/*.mod *.out *.file analyze *~ 2>/dev/null

Files="geometry.f90 grid_generation/matmul.f90 grid_generation/simple_rotations.f90 grid_generation/spherical.f90"
Files=$Files" grid_generation/projections.f90 grid_generation/matrix_rotation.f90 grid_generation/grid_generator.f90"
Files=$Files" parallel_cubic.f90 metrics.f90 grid_var.f90 derivatives.f90 interpolation.f90 func_var.f90 messenger.f90 diagnostic.f90 printer.f90 schemes.f90"
Files=$Files" Solver_shallow_water.f90"

netcdf="/data4t/avolkov/util/netcdf-2016Jan-13.1"
netcdf="/home/sasha/netcdf"

if [[ $grid_type == 1 ]]; then
	grid="equiang"
elif [ $grid_type == 0 ]&&[ $rescale == 0 ]; then
	grid="simple"
elif [ $grid_type == 0 ] && [ $rescale == 1 ]; then
	grid="tan"
fi

DIRECTORY="datFiles/$(( 2*$dim ))/$grid"
PIC="datFiles/$(( 2*$dim ))/pic"

if [ ! -d "$DIRECTORY" ]; then
	mkdir $DIRECTORY
fi
if [ ! -d "$PIC" ]; then
	mkdir $PIC
fi
if [ -d "$DIRECTORY" ]; then
	echo $DIRECTORY
fi

 # -check all -traceback -ftrapuv
mpiifort -openmp -O3 $Files -module mod_files -I grid_generation -I $netcdf/inc -L $netcdf/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm 2> err.file
# /home/sasha/Fortran/Comands/./compo geometry.o conformal.o matmul.o morphism.o grid_generator.o data_analyzer.o spherical.o main.o
	CompStatus=$?
	echo "compilation status" $CompStatus

if [[ `grep -c error err.file` > 0 ]]; then
	echo "Look for" `grep -c error err.file` "errors in err.file"
	echo `grep -c warning err.file` "warnings"
else

	echo `grep -c error err.file` "errors"
	echo `grep -c warning err.file` "warnings"

	if [[ $1 != "compile" ]]; then
		export OMP_NUM_THREADS=4
		mpirun -n $1 ./a.out

		echo "Regridding"
		cd regrid
		./run_intel_reg.sh

		echo "Plotting"
		cd ..
		./plotscript.sh
	fi

fi
