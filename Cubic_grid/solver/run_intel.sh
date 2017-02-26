#/bin/bash

line=$(cat init | sed 's/ //g')
args=(${line//,/ })
dim="${args[0]}"
rescale="${args[4]}"
grid_type="${args[5]}"
step="${args[6]}"

if [ ! -d mod_files ]; then
	mkdir mod_files
fi

rm -rf *.o mod_files/*.mod *.out *.file analyze *~ 2>/dev/null

gg=grid_generation

Files="$gg/matmul.f90 $gg/simple_rotations.f90 $gg/spherical.f90"
Files=$Files" $gg/projections.f90 $gg/matrix_rotation.f90 $gg/grid_generator.f90"
Files=$Files" geometry.f90 parallel_cubic.f90 grid_var.f90 metrics.f90 derivatives.f90 interpolation.f90 func_var.f90"
Files=$Files" messenger.f90 diagnostic.f90 printer.f90 methods.f90 Solver_shallow_water.f90"

netcdf="/home/sasha/netcdf"
if [ ! -d "$netcdf" ]; then
	netcdf="/data4t/avolkov/util/netcdf-2016Jan-13.1"
fi

if [[ $grid_type == 1 ]]; then
	grid="equiang"
elif [ $grid_type == 0 ]&&[ $rescale == 0 ]; then
	grid="simple"
elif [ $grid_type == 0 ] && [ $rescale == 1 ]; then
	grid="tan"
fi

DIRECTORY="datFiles/$(( 2*$dim ))"
SubDIRECTORY="datFiles/$(( 2*$dim ))/$(( 2*$step ))th"
SSubDIRECTORY="datFiles/$(( 2*$dim ))/$(( 2*$step ))th/$grid"
PIC="datFiles/pic/$(( 2*$step ))th"

if [ ! -d datFiles ]; then
	mkdir datFiles
fi
if [ ! -d "$DIRECTORY" ]; then
	mkdir $DIRECTORY
fi
if [ ! -d "$SubDIRECTORY" ]; then
	mkdir $SubDIRECTORY
fi
if [ ! -d "$SSubDIRECTORY" ]; then
	mkdir $SSubDIRECTORY
fi
if [ ! -d "$PIC" ]; then
	mkdir $PIC
fi
if [ -d "$SubDIRECTORY" ]; then
	echo $SSubDIRECTORY
fi

if [[ $1 == "compile" ]]; then
mpiifort -check all -traceback -ftrapuv -openmp $Files -module mod_files -I $gg -I $netcdf/inc -L $netcdf/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm 2> err.file
elif [[ $1 != "compile" ]] ; then
mpiifort -openmp -O3 $Files -module mod_files -I $gg -I $netcdf/inc -L $netcdf/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm 2> err.file
# mpiifort -check all -traceback -ftrapuv -g $Files -module mod_files -I $gg -I $netcdf/inc -L $netcdf/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -lm 2> err.file
fi
# /home/sasha/Fortran/Comands/./compo geometry.o conformal.o matmul.o morphism.o grid_generator.o data_analyzer.o spherical.o main.o
	CompStatus=$?
	echo "compilation status" $CompStatus

if [[ `grep -c error err.file` > 0 ]]; then
	echo "Look for" `grep -c error err.file` "errors in err.file"
	echo `grep -c warning err.file` "warnings"
	echo `grep -c ifort err.file` "ifort"
else

	echo `grep -c error err.file` "errors"
	echo `grep -c warning err.file` "warnings"
	echo `grep -c ifort err.file` "ifort"

	if [[ $1 != "compile" ]]; then
		export OMP_NUM_THREADS=$2
		mpirun -n $1 ./a.out

		echo "Regridding"
		cd regrid
		./run_intel_reg.sh

	fi

fi
