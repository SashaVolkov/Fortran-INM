program solver

	use sphere_geometry, Only: geometry
	use grid_var, Only: g_var
	use metrics, Only: metric
	use parallel_cubic, Only: parallel
	use func_var, Only: f_var
	use printer_ncdf, Only: printer
	use methods, Only: method
	use diagnostic_mod, Only: diagnostic
	use messenger, Only: message
	use interpolation, Only: interp
	use omp_lib
	use mpi

implicit none


!variables
	real(8) r_sphere, g, pi, step, omega_cor, height, dt, start_init, end_init
	integer(4) dim, space_step, Tmax, time, speedup, Wid, grid_id, xid, yid, faceid, ncid, ncid_gr, rescale, face, grid_type, x, y, flag

	integer(4) status(MPI_STATUS_SIZE), ier, id, np, numthreads

	Type(geometry) :: geom
	Type(f_var) :: var, var_prev
	Type(g_var) :: grid
	Type(metric) :: metr
	Type(parallel) :: paral
	Type(printer) :: printer_nc
	Type(method) :: met
	Type(diagnostic) :: diagn
	Type(message) :: msg
	Type(interp) :: inter


!definition
	r_sphere= 6371220d0;  g = 9.80616
	pi = 314159265358979323846d-20;  omega_cor = 7292d-2
	height = 100.0;  dt = 30.0;  flag = 0
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

	call paral.init(dim, space_step+1, np, id) ! Yep, that's right step+1 is correct, interpolation and transformation matrix need more than step points
	call geom.init(r_sphere, pi)
	call metr.init(paral)
	call grid.init(geom, paral, metr, omega_cor, g, dt, rescale, grid_type)
	call var.init(paral, metr, height)
	call var_prev.init(paral, metr, height)
	call var_prev.start_conditions(metr, geom)
	call diagn.init( grid, paral, Tmax, id)
	call met.init(var_prev, grid, space_step)
	call msg.init(grid_type)
	call inter.init(grid, 4)

	call diagn.L_norm(var_prev, grid, time)

	call printer_nc.init(dim, Tmax, speedup, time, grid_id, ncid, ncid_gr, rescale, grid_type)
	call printer_nc.to_print(var_prev, diagn, 0, speedup, ncid, id)
	if(id == 0) call printer_nc.print_grid(grid, grid_id, ncid_gr)


	do time = 1, Tmax
		! call met.Euler(var, var_prev, metr, inter, paral, msg)
		! call met.Predictor_corrector(var, var_prev, metr, inter, paral, msg)
		call met.RungeKutta(var, var_prev, metr, inter, paral, msg)
			if(mod(time, speedup) == 0) then
				call diagn.Courant(var_prev, metr, time)
				call printer_nc.to_print(var_prev, diagn, time, speedup, ncid, id)
				call diagn.L_norm(var_prev, grid, time)
			end if
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
		print '(" max/min = ", f10.2)', grid.max_to_min
		print '(" np = ", I5)', np
		print '(" time = ", f10.2, " sec")', end_init - start_init
	end if


	call grid.deinit()
	call var.deinit()
	call var_prev.deinit()
	call printer_nc.deinit(ncid, ncid_gr)
	call diagn.deinit()
	call met.deinit()
	call metr.deinit()


	call MPI_FINALIZE(ier)

end program