program solver

	use grid_generator_solver, Only: grid

!variables
	real(8) t(2), r_sphere
	integer(4) dimention
	real(8), Allocatable :: r_out(:, :, :, :)
	Type(grid) :: generator
!definition

	r_sphere = 1d0
	dimention = 400
	Allocate(r_out(1:6, -dimention:dimention, -dimention:dimention, 1:3))

!subroutines calls
	call cpu_time(t(1))
	call generator.conformal_cubed_sphere(dimention,dimention, r_sphere, r_out)
	call cpu_time(t(2))
	print '("Time of generation = ", f6.3, " sec")', t(2) - t(1)









end program