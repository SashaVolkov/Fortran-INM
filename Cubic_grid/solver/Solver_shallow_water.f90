program solver

	use grid_generator_solver, Only: grid
	use special_variables, Only: variables

implicit none

!variables
	real(8) r_sphere, g, pi, step, omega_cor
	integer(4) dim ! dimension
	real(8), Allocatable :: grid_points(:, :, :, :)
	Type(grid) :: generator
	Type(variables) :: var, var_prev


!definition
	r_sphere= 6371220d0; g = 980616d-5; pi = 314159265358979323846d-20; omega_cor = 7292d-2
	dim = 400
	step = 2*pi*r_sphere/(8d0*dim)
	Allocate(grid_points(1:6, -dim:dim, -dim:dim, 1:2)) ! face_id, 2dim*2dim, latitude(theta) & longitude(lambda)


!subroutines calls
	call generator.conformal_cubed_sphere(dim,dim, r_sphere, grid_points)

	call var.init(grid_points, dim, omega_cor, r_sphere, g)
	call var_prev.init(grid_points, dim, omega_cor, r_sphere, g)
	call var_prev.start_conditions(grid_points, dim)



	print '(" Theta = ", f7.4, " Phi = ", f7.4)', grid_points(6, 0, 0, :)
	print '(" Theta = ", f7.4, " Phi = ", f7.4)', grid_points(4, 0, 0, :)
	print '(" Theta = ", f7.4, " Phi = ", f7.4)', grid_points(2, 0, 0, :)
	print '(" Theta = ", f7.4, " Phi = ", f7.4)', grid_points(1, 0, 0, :)
	print '(" Grid step = ", f10.2, " m")', step





end program