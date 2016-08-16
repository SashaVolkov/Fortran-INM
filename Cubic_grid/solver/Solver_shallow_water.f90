program solver

	use sphere_geometry, Only: geometry
	use grid_var, Only: g_var
	use parallel_cubic, Only: parallel
	use func_var, Only: f_var
	use printer_ncdf, Only: printer
	use schemes, Only: schema
	use diagnostic_mod, Only: diagnostic
	use messenger, Only: message
	use interpolation, Only: interp
! 	use omp_lib
	use mpi

implicit none


!variables
	real(8) r_sphere, g, pi, step, omega_cor, height, dt, start_init, end_init
	integer(4) dim, gr_step, Tmax, time, speedup, Wid, xid, yid, faceid, ncid, rescale, face

	integer(4) status(MPI_STATUS_SIZE), ier, id, np, numthreads

	Type(geometry) :: geom
	Type(f_var) :: var, var_prev
	Type(g_var) :: grid
	Type(parallel) :: paral
	Type(printer) :: printer_nc
	Type(schema) :: sch
	Type(diagnostic) :: diagn
	Type(message) :: msg
	Type(interp) :: inter


!definition
	r_sphere= 6371220d0;  g = 980616d-5
	pi = 314159265358979323846d-20;  omega_cor = 7292d-2
	dim = 25;  gr_step = 2;  height = 100.0
	step = 2*pi*r_sphere/(8d0*dim)

	Tmax = 250;  speedup = 250;  dt = 5d0
	rescale = 0 ! 0-simple, 1-tan, 2-pow(4/3)q
!480000

	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)

!subroutines calls

	start_init = MPI_Wtime()

	call paral.init(dim, gr_step, np, id)
	call geom.init(r_sphere, pi)
	call grid.init(geom, paral, omega_cor, g, dt, rescale)
	call var.init(paral, height)
	call var_prev.init(paral, height)
	call sch.init(var_prev, grid)
	call msg.init()
	call inter.init(grid, 1)

	call var_prev.start_conditions()

	call printer_nc.init(dim, Tmax, speedup, time, Wid, xid, yid, faceid, ncid, rescale)
	call printer_nc.to_print(var_prev, 0, speedup, Wid, ncid, id)
	call diagn.init( grid, paral, Tmax, rescale, id)


	! do time = 1, Tmax
	! 	call sch.RungeKutta(var, var_prev, grid)
	! 	call msg.msg(var_prev, paral)
	! 	call diagn.L_norm(var_prev, grid, time)
	! 	call diagn.Courant(var_prev, grid, time)
	! 		if(mod(time, speedup) == 0) call printer_nc.to_print(var_prev, time, speedup, Wid, ncid, id)
	! 		if(mod(time, Tmax/10) == 0 .and. id == 0) then
	! 			end_init = MPI_Wtime()
	! 			print '(I3, "% Done time = ", f7.2, " sec")', time*100/Tmax, end_init - start_init
	! 		end if
	! end do

	end_init = MPI_Wtime()

	if(id == 0) then
		print '(" Grid step =  ", f10.2, " m")', step
		print '(" Y max step = ", f10.2, " m")', grid.dy_max * r_sphere
		print '(" Y min/max = ", f6.4)', grid.dy_min/grid.dy_max
		print '(" X max/min = ", f6.4)', grid.dx_max/grid.dx_min
		! print '(" latlon = ", f8.3, f8.3)', grid.latlon(1, 1, 1, 2) * 180.0/pi
		! print '(" latlon = ", f8.3, f8.3)', grid.latlon(1, 0, 1, 2) * 180.0/pi
		print '(" np = ", I5)', np
		print '(" time = ", f10.2, " sec")', end_init - start_init
	end if

	! print *, "threads = ", omp_get_num_threads()

	call grid.deinit()
	call var.deinit()
	call var_prev.deinit()
	call printer_nc.deinit()
	call diagn.deinit()
	call sch.deinit()


	call MPI_FINALIZE(ier)

end program