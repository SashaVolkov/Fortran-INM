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
	Real(8) r_sphere, g, pi, step, omega_cor, height, dt, start_init, end_init
	Integer(4) dim, space_step, Tmax, time, speedup, rescale, face, grid_type, x, y, flag, test, exit_flag, i
	Real(8), Allocatable :: cycle_time(:)

	Integer(4) status(MPI_STATUS_SIZE), ier, id, np, numthreads

	Type(geometry) :: geom
	Type(f_var) :: var, var_prev
	Type(g_var) :: grid
	Type(metric) :: metr
	Type(parallel) :: paral
	Type(printer) :: printer_nc
	Type(method) :: meth
	Type(diagnostic) :: diagn
	Type(message) :: msg
	Type(interp) :: inter


!definition
	r_sphere= 6371220d0
	pi = 314159265358979323846d-20;  omega_cor = 7292d-8 !13d-5 ! 
	height = 100d0;  flag = 0
! r_sphere= 1d0
	! rescale  0-simple, 1-tan, 2-pow(4/3)q
	! grid_type  0 - conformal, 1 - equiangular


	call MPI_Init(ier)
	call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
	call MPI_Comm_size(MPI_COMM_WORLD,np,ier)

!subroutines calls

	open(9,file='init')
		read(9, *) dim, Tmax, dt, speedup, rescale, grid_type, space_step, test
	close(9)

	if ( id == 0 ) Allocate(cycle_time(1:Tmax))


	call paral.init(dim, 2*space_step, np, id) ! Yep, that's right step+1 is correct, interpolation and transformation matrix need more than step points
	call msg.init(grid_type, paral)
	call geom.init(r_sphere, pi)
	call grid.init(geom, paral, rescale, grid_type)
	call metr.init(grid)
	call var.init(metr, test)
	call var_prev.init(metr, test)
	call var_prev.start_conditions(metr, geom, omega_cor)
	call diagn.init(var_prev, grid, Tmax, id, dt)
	call meth.init(var_prev, space_step, Tmax, dt)
	call inter.init(metr, 4)
	call printer_nc.init(grid, Tmax, speedup, time, rescale, grid_type, test)


	call printer_nc.to_print(var_prev, diagn, 0)
	call printer_nc.print_grid(grid)

	call MPI_Barrier(MPI_COMM_WORLD, ier)
	start_init = MPI_Wtime();  cycle_time = 0d0

	do time = 1, Tmax
! 		call meth.Euler(var, var_prev, metr, inter, msg)
! 		call meth.Predictor_corrector(var, var_prev, metr, inter, msg)
		if ( id == 0 ) cycle_time(time) = MPI_Wtime()
		call meth.RungeKutta(var, var_prev, metr, inter, msg, time)

		do i = 0, np-1
			if ( isnan(var.h_height(var.ns_x,var.ns_y,1)) ) exit_flag = 1
			call MPI_Bcast(exit_flag, 1, MPI_INT, i, MPI_COMM_WORLD)
		end do
		if ( exit_flag == 1  .and. id == 0 ) print '(I3, "% Done time = ", f7.2, " sec")', time*100/Tmax, end_init - start_init
		if ( exit_flag == 1 ) exit

		if ( id == 0 ) cycle_time(time) = MPI_Wtime() - cycle_time(time)
				call diagn.Courant(metr, var_prev, time)
			if(mod(time, speedup) == 0) then
				call printer_nc.to_print(var_prev, diagn, time)
			end if
		end_init = MPI_Wtime()
			if(mod(10*time, Tmax) == 0 ) call MPI_Barrier(MPI_COMM_WORLD, ier)
			if(mod(10*time, Tmax) == 0 .and. id == 0) then
				print '(I3, "% Done time = ", f7.2, " sec")', time*100/Tmax, end_init - start_init
			end if
	end do


	end_init = MPI_Wtime()

	step = 2*pi*r_sphere/(8d0*dim)
	if(id == 0) then
		print '(" Grid step =  ", f10.2, " m")', step
		print '(" Grid step =  ", f10.2, " m")', grid.delta_on_cube
		print '(" max/min = ", f10.2)', grid.max_to_min
		print '(" min = ", f10.2)', grid.min
		print '(" np = ", I5)', np
		print '(" threads = ", I5)', OMP_GET_NUM_THREADS()
		print '(" max = ", f10.6, " sec")', maxval(cycle_time)
		print '(" min = ", f10.6, " sec")', minval(cycle_time)
		print '(" average = ", f10.6, " sec")', sum(cycle_time)/time
		print '(" full = ", f10.2, " sec")', sum(cycle_time)
		print '(" msg time = ", f10.6, " sec")', sum(meth.msg_time)
		print '(" msg time average = ", f10.8, " sec")', sum(meth.msg_time)/(5*Tmax)
		print '(" time = ", f10.2, " sec")', end_init - start_init
	end if


	call grid.deinit()
	call var.deinit()
	call var_prev.deinit()
	call printer_nc.deinit()
	call diagn.deinit()
	call meth.deinit()
	call metr.deinit()


	call MPI_FINALIZE(ier)

end program