program solver

	use grid_var, Only: g_var
	use func_var, Only: f_var
	use printer_ncdf, Only: printer
	use schemes, Only: schema
	use diagnostic_mod, Only: diagnostic

implicit none

!variables
	real(8) r_sphere, g, pi, step, omega_cor, height, dt
	integer(4) dim, gr_step, Tmax, time, speedup, Wid, xid, yid, ncid(1:6), rescale

	Type(f_var) :: var, var_prev
	Type(g_var) :: grid
	Type(printer) :: printer_nc
	Type(schema) :: sch
	Type(diagnostic) :: diagn


!definition
	r_sphere= 6371220d0;  g = 980616d-5
	pi = 314159265358979323846d-20;  omega_cor = 7292d-2
	dim = 25;  gr_step = 1;  height = 100.0
	step = 2*pi*r_sphere/(8d0*dim)

	Tmax = 80000;  speedup = 80;  dt = 500d0
	rescale = 1 ! 0-simple, 1-tan



!subroutines calls

	print '(" init")'


	call grid.init(dim, gr_step, omega_cor, r_sphere, g, dt, rescale)
	dim = grid.dim
	call var.init(dim, gr_step, height)
	call var_prev.init(dim, gr_step, height)
	call diagn.init(grid, Tmax, rescale, pi)

	call var_prev.start_conditions()

	call printer_nc.init(dim, Tmax, speedup, time, Wid, xid, yid, ncid)
	call printer_nc.to_print(var_prev, dim, 0, speedup, Wid, ncid)


	do time = 1, Tmax
		call sch.Linear(var, var_prev, grid)
		if(mod(time, speedup) == 0) call diagn.CFL(var_prev, grid, time)
		if(mod(time, speedup) == 0) call printer_nc.to_print(var_prev, grid.dim, time, speedup, Wid, ncid)
	end do


	print '(" Grid step = ", f10.2, " m")', step
	print '(" X min step = ", f10.2, " m")', grid.dx_min
! 	print '(" Y min step = ", f10.2, " m")', grid.dy_min
! 	print '(" X max step = ", f10.2, " m")', grid.dx_max
	print '(" Y max step = ", f10.2, " m")', grid.dy_max
	print '(" Y max/min = ", f10.2)', grid.dy_max/grid.dy_min
	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,dim,dim)
	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,dim,0)
	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,0,0)

! 	! print '(" Coriolis = ", f10.2)', var.f_cor(-90, -100, 1)


	call grid.deinit()
	call var.deinit()
	call var_prev.deinit()
	call printer_nc.deinit()
	call diagn.deinit()


end program