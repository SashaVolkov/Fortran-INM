program solver

	use sphere_geometry, Only: geometry
	use grid_var, Only: g_var
	use parallel_cubic, Only: parallel
	use func_var, Only: f_var
	use printer_ncdf, Only: printer
	use schemes, Only: schema
	use diagnostic_mod, Only: diagnostic

implicit none

	include"mpif.h"

!variables
	real(8) r_sphere, g, pi, step, omega_cor, height, dt
	integer(4) dim, gr_step, Tmax, time, speedup, Wid, xid, yid, ncid(1:6), rescale, face

	integer(4) status(MPI_STATUS_SIZE), ier, id, np

	Type(geometry) :: geom
	Type(f_var) :: var(1:6), var_prev(1:6)
	Type(g_var) :: grid
	Type(parallel) :: par
	Type(printer) :: printer_nc
	Type(schema) :: sch
	Type(diagnostic) :: diagn


!definition
	r_sphere= 6371220d0;  g = 980616d-5
	pi = 314159265358979323846d-20;  omega_cor = 7292d-2
	dim = 400;  gr_step = 1;  height = 100.0
	step = 2*pi*r_sphere/(8d0*dim)

	Tmax = 6000;  speedup = 40;  dt = 400d0
	rescale = 1 ! 0-simple, 1-tan, 2-pow(4/3)


	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)

!subroutines calls


	call geom.init(r_sphere, pi)
	call grid.init(geom, dim, gr_step, omega_cor, g, dt, rescale)
	dim = grid.dim

	call par.init(dim, gr_step, np, id)

! 	do face = 1, 6
! 		call var(face).init(par, dim, gr_step, height, face)
! 		call var_prev(face).init(par, dim, gr_step, height, face)
! 		call var_prev(face).start_conditions()
! 	end do

! ! 	print *, var(1).ns_x, var(1).ns_y, var(1).nf_x, var(1).nf_y
! ! 	print *, var_prev(1).ns_x, var_prev(1).ns_y, var_prev(1).nf_x, var_prev(1).nf_y
! ! 	! print *, par.ns_xy, par.nf_xy

! 	call printer_nc.init(dim, Tmax, speedup, time, Wid, xid, yid, ncid, rescale)
! 	call printer_nc.to_print(var_prev, 0, speedup, Wid, ncid, id)
! ! 	diagn.init( grid, Tmax, rescale)

! ! 			print '(" calc")'

! 			do time = 1, Tmax
! 				call sch.Linear(var, var_prev, grid)
! ! ! 				if(mod(time, speedup) == 0) call diagn.Courant(var_prev, grid, time)
! ! 				if(mod(time, speedup) == 0) call diagn.L_norm(var_prev, grid, time)
! 				if(mod(time, speedup) == 0) call printer_nc.to_print(var_prev, time, speedup, Wid, ncid, id)
! 				! call MPI_Barrier(MPI_COMM_WORLD, ier)
! 			end do

! 			print *, np
			print '(" Grid step = ", f10.2, " m")', step
! 			print '(" X min step = ", f10.2, " m")', grid.dx_min
! 			print '(" Y max step = ", f10.2, " m")', grid.dy_max
			print '(" Y max/min = ", f10.2)', grid.dy_max/grid.dy_min

! 			print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,dim,dim)
! 			print '(" ")'
! 			print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,12,0)
! 			print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,0,12)
! 			print '(" ")'
! 			print '(" step = ", f10.2, f10.2, f10.2, f10.2)', grid.h_dist(:,0,0)



	call grid.deinit()
	do face = 1, 6
		call var(face).deinit()
		call var_prev(face).deinit()
	end do
	call printer_nc.deinit()
	call diagn.deinit()


	call MPI_FINALIZE(ier)

end program