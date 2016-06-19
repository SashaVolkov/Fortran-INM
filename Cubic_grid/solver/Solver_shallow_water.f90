program solver

	use sphere_geometry, Only: geometry
	use grid_var, Only: g_var
	use parallel_cubic, Only: parallel
	use func_var, Only: f_var
	use printer_ncdf, Only: printer
	use schemes, Only: schema
	use diagnostic_mod, Only: diagnostic
	use messenger, Only: message

implicit none

	include"mpif.h"

!variables
	real(8) r_sphere, g, pi, step, omega_cor, height, dt
	integer(4) dim, gr_step, Tmax, time, speedup, Wid, xid, yid, ncid(1:6), rescale, face

	integer(4) status(MPI_STATUS_SIZE), ier, id, np

	Type(geometry) :: geom
	Type(f_var) :: var, var_prev
	Type(g_var) :: grid
	Type(parallel) :: paral
	Type(printer) :: printer_nc
	Type(schema) :: sch
	Type(diagnostic) :: diagn
	Type(message) :: msg


!definition
	r_sphere= 6371220d0;  g = 980616d-5
	pi = 314159265358979323846d-20;  omega_cor = 7292d-2
	dim = 32;  gr_step = 2;  height = 100.0
	step = 2*pi*r_sphere/(8d0*dim)

	Tmax = 8000;  speedup = 40;  dt = 400d0
	rescale = 1 ! 0-simple, 1-tan, 2-pow(4/3)


	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)

!subroutines calls


	call geom.init(r_sphere, pi)
	call grid.init(geom, dim, gr_step, omega_cor, g, dt, rescale)
	dim = grid.dim
	call paral.init(dim, gr_step, np, id)
	call var.init(paral, gr_step, height)
	call var_prev.init(paral, gr_step, height)
	call msg.init()

	call var_prev.start_conditions()

	call printer_nc.init(dim, Tmax, speedup, time, Wid, xid, yid, ncid, rescale)
	call printer_nc.to_print(var_prev, 0, speedup, Wid, ncid, id)
! 	diagn.init( grid, Tmax, rescale)


	do time = 1, 1
		call sch.Linear(var, var_prev, grid)
! 				if(mod(time, speedup) == 0) call diagn.Courant(var_prev, grid, time)
! 				if(mod(time, speedup) == 0) call diagn.L_norm(var_prev, grid, time)
		call msg.msg(var_prev, paral)
		if(mod(time, speedup) == 0) call printer_nc.to_print(var_prev, time, speedup, Wid, ncid, id)
	end do


	if(id == 0) then
		print '(" Grid step =  ", f10.2, " m")', step
		print '(" Y max step = ", f10.2, " m")', grid.dy_max
		print '(" Y max/min = ", f5.2)', grid.dy_max/grid.dy_min
		print '(" X max/min = ", f5.2)', grid.dx_max/grid.dx_min
		! print '(" latlon = ", f10.2, f10.2)', grid.points_latlon(:, dim+1, 1, 4) * 180.0/pi
		! print '(" latlon = ", f10.2, f10.2)', grid.points_latlon(:, dim+1, 2*dim+1, 4) * 180.0/pi
		print '(" np = ", I7)', np
	end if



	call grid.deinit()
	call var.deinit()
	call var_prev.deinit()
	call printer_nc.deinit()
	call diagn.deinit()


	call MPI_FINALIZE(ier)

end program