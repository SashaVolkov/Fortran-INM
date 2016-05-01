program solver

	use grid_var, Only: g_var
	use func_var, Only: f_var
	use printer_ncdf, Only: printer
	use schemes, Only: schema
	use diagnostic_mod, Only: diagnostic

implicit none

!variables
	real(8) r_sphere, g, pi, step, omega_cor, height, dt
	integer(4) dim, gr_step, Tmax, time, speedup, Wid, xid, yid, ncid(1:6)
	! real(8), Allocatable :: grid_points(:, :, :, :)

	Type(f_var) :: var, var_prev
	Type(g_var) :: grid
	Type(printer) :: printer_nc
	Type(schema) :: sch
	Type(diagnostic) :: diagn


!definition
	r_sphere= 6371220d0;  g = 980616d-5
	pi = 314159265358979323846d-20;  omega_cor = 7292d-2
	dim = 400;  gr_step = 1;  height = 100.0
	step = 2*pi*r_sphere/(8d0*dim)

	Tmax = 40000;  speedup = 80;  dt = 500d0

	! Allocate(grid_points(1:2, -dim:dim, -dim:dim, 1:6)) ! latitude & longitude, 2dim*2dim, face_id


!subroutines calls

	print '(" init")'


	call grid.init(dim, gr_step, omega_cor, r_sphere, g, dt)
! 	call var.init(dim, gr_step, height)
! 	call var_prev.init(dim, gr_step, height)
! 	call var_prev.start_conditions()

! 	call printer_nc.init(dim, Tmax, speedup, time, Wid, xid, yid, ncid)
! 	call printer_nc.to_print(var_prev, dim, 0, speedup, Wid, ncid)

! 	call diagn.init(grid, Tmax)

	! Main cycle
! 	do time = 1, Tmax
! 		call sch.Linear(var, var_prev, grid)
! 		if(mod(time, speedup) == 0) call diagn.CFL(var_prev, grid, time)
! 		if(mod(time, speedup) == 0) call printer_nc.to_print(var_prev, dim, time, speedup, Wid, ncid)
! 	end do
	! End of cycle

	print '(" Grid step = ", f10.2, " m")', step
	print '(" X min step = ", f10.2, " m")', grid.dx_min
	print '(" Y min step = ", f10.2, " m")', grid.dy_min
	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,dim,dim)
	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,dim,dim-1)
	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,dim,dim-2)
	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,dim,dim-3)
	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,dim,0)
! 	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,-dim,-dim)
! 	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,-300,-300)
! 	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,-100,-100)
! 	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,12,12)
	print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,0,0)
	! print '(" Coriolis = ", f10.2)', var.f_cor(-90, -100, 1)


	call grid.deinit()
	call var.deinit()
	call var_prev.deinit()
	call printer_nc.deinit()
	call diagn.deinit()


end program