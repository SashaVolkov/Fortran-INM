program solver

	use grid_generator_solver, Only: grid
	use special_variables, Only: variables
	use printer_ncdf, Only: printer
	use schemes, Only: schema

implicit none

!variables
	real(8) r_sphere, g, pi, step, omega_cor, height
	integer(4) dim, gr_step ! dimension
	real(8), Allocatable :: grid_points(:, :, :, :)
	Type(grid) :: generator
	Type(variables) :: var, var_prev
	Type(printer) :: printer_nc
	Type(schema) :: sch


!definition
	r_sphere= 6371220d0;  g = 980616d-5
	pi = 314159265358979323846d-20;  omega_cor = 7292d-2
	dim = 400; gr_step = 1;  height = 100.0
	step = 2*pi*r_sphere/(8d0*dim)

	Allocate(grid_points(1:2, -dim:dim, -dim:dim, 1:6)) ! latitude & longitude, 2dim*2dim, face_id


!subroutines calls
	call generator.conformal_cubed_sphere(dim,dim, r_sphere, grid_points)

	call var.init(grid_points, dim, gr_step, omega_cor, r_sphere, g, height)
	call var_prev.init(grid_points, dim, gr_step, omega_cor, r_sphere, g, height)
	call var_prev.start_conditions(dim)

	call printer_nc.init()
	call printer_nc.to_print()

	! Main cycle
	call sch.Linear(var, var_prev)
	call printer_nc.to_print()
	! End of cycle

	print '(" Grid step = ", f10.2, " m")', step


	call var.deinit()
	call var_prev.deinit()
	call printer_nc.deinit()

end program