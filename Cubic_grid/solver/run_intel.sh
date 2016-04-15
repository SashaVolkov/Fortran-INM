#/bin/bash

rm -rf *.o *.mod *.out *.file analyze grid *~ 2>/dev/null



mpiifort geometry.f90 matmul.f90 simple_rotations.f90 spherical.f90 projections.f90 matrix_rotation.f90 grid_generator.f90 Solver_shallow_water.f90 2> err.file
# /home/sasha/Fortran/Comands/./compo geometry.o conformal.o matmul.o morphism.o grid_generator.o data_analyzer.o spherical.o main.o
	echo "compilation status" $?

if [[ `grep -c error err.file` > 0 ]]; then
	echo "Look for" `grep -c error err.file` "errors in err.file"
	echo `grep -c warning err.file` "warnings"
else

	echo `grep -c error err.file` "errors"
	echo `grep -c warning err.file` "warnings"

	time mpiexec -n 1 ./a.out

fi


