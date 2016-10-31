program solver

	use sphere_geometry, Only: geometry
	use grid_var, Only: g_var
	use metrics, Only: metric
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
	integer(4) dim, space_step, Tmax, time, speedup, Wid, grid_id, xid, yid, faceid, ncid, ncid_gr, rescale, face, grid_type

	integer(4) status(MPI_STATUS_SIZE), ier, id, np, numthreads

	Type(geometry) :: geom
	Type(f_var) :: var, var_prev
	Type(g_var) :: grid
	Type(metric) :: metr
	Type(parallel) :: paral
	Type(printer) :: printer_nc
	Type(schema) :: sch
	Type(diagnostic) :: diagn
	Type(message) :: msg
	Type(interp) :: inter


!definition
	r_sphere= 6371220d0;  g = 9.80616
	pi = 314159265358979323846d-20;  omega_cor = 7292d-2
	height = 100.0;  dt = 30.0
! r_sphere= 1d0
	! rescale  0-simple, 1-tan, 2-pow(4/3)q
	! grid_type  0 - conformal, 1 - equiangular


	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)

!subroutines calls

	open(9,file='init')
		read(9, *) dim, Tmax, dt, speedup, rescale, grid_type, space_step
	close(9)

	start_init = MPI_Wtime()

	call paral.init(dim, space_step, np, id)
	call geom.init(r_sphere, pi)
	call metr.init(paral)
	call grid.init(geom, paral, metr, omega_cor, g, dt, rescale, grid_type)
	call var.init(paral, metr, height)
	call var_prev.init(paral, metr, height)
	call var_prev.start_conditions(metr, geom)
	call diagn.init( grid, paral, Tmax, id)
	call sch.init(var_prev, grid, space_step)
	call msg.init(grid_type)
	call inter.init(grid, 2)


	call printer_nc.init(dim, Tmax, speedup, time, Wid, grid_id, ncid, ncid_gr, rescale, grid_type)
	call printer_nc.to_print(var_prev, 0, speedup, Wid, ncid, id)
	if(id == 0) call printer_nc.print_grid(grid, grid_id, ncid_gr)


	do time = 1, Tmax
		! call sch.Linear(var, var_prev, grid, metr, inter, paral, msg)
! 		call sch.RungeKutta(var, var_prev, grid, metr, inter, paral, msg)
		call sch.INM_sch(var, var_prev, grid, metr, inter, paral, msg)
		call diagn.L_norm(var_prev, grid, time)
		call diagn.Courant(var_prev, grid, metr, time)
			if(mod(time, speedup) == 0) call printer_nc.to_print(var_prev, time, speedup, Wid, ncid, id)
			if(mod(10*time, Tmax) == 0 .and. id == 0) then
				end_init = MPI_Wtime()
				print '(I3, "% Done time = ", f7.2, " sec")', time*100/Tmax, end_init - start_init
			end if
	end do


	end_init = MPI_Wtime()

	step = 2*pi*r_sphere/(8d0*dim)
	if(id == 0) then
		print '(" Grid step =  ", f10.2, " m")', step
		print '(" Grid step =  ", f10.2, " m")', grid.delta_on_cube
! 		print '(" Y max step = ", f10.2, " m")', grid.dy_max
! 		print '(" Y min step = ", f10.2, " m")', grid.dy_min
! 		print '(" Y min/max = ", f6.4)', grid.dy_min/grid.dy_max
! 		print '(" X max/min = ", f6.4)', grid.dx_max/grid.dx_min
		! print '(" latlon = ", f8.3, f8.3)', grid.latlon(1, 1, 1, 2) * 180.0/pi
		! print '(" latlon = ", f8.3, f8.3)', grid.latlon(1, 0, 1, 2) * 180.0/pi
		print '(" np = ", I5)', np
		print '(" time = ", f10.2, " sec")', end_init - start_init
	end if


	call grid.deinit()
	call var.deinit()
	call var_prev.deinit()
	call printer_nc.deinit(ncid, ncid_gr)
	call diagn.deinit()
	call sch.deinit()
	call metr.deinit()


	call MPI_FINALIZE(ier)

end program