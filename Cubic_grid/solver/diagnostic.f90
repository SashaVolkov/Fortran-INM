module diagnostic_mod

	use func_var, Only: f_var
	use grid_var, Only: g_var
	use metrics, Only: metric
	use mpi
	use omp_lib


	implicit none
	Private
	Public :: diagnostic


	Type diagnostic

		Real(8), Allocatable :: CFL(:,:,:)
		Real(4), Allocatable :: square(:, :)
		Integer(4) Tmax, dim, step, flag
		Real(8) convert_time, L10, L20, L_inf0, dh, dt, center(1:2)


		CONTAINS
			Procedure, Public :: init => init
			Procedure, Private :: alloc => alloc
			Procedure, Public :: deinit => deinit
			Procedure, Public :: Courant => Courant
			Procedure, Public :: histogram => hist_generation

	End Type


CONTAINS



	Subroutine init(this, func, grid, Tmax, id, dt)

		Class(diagnostic) :: this
		Class(f_var) :: func
		Class(g_var) :: grid

		Integer(4), intent(in) :: Tmax, id
		Real(8), intent(in) :: dt
		character(32) istring, istring1

		this.Tmax = Tmax;  this.dim = func.dim;  this.step = func.step-1
		this.dt = dt;  this.dh = func.delta_on_cube
		this.convert_time = this.dt/(60d0*60d0*24d0)
		this.flag = 1

		write(istring, *) 2*this.dim
		write(istring1, *) 2*this.step

		istring = trim(adjustl(istring))//'/'//trim(adjustl(istring1))//'th'

		call this.alloc(func)

		if (func.grid_type == 1) then
			istring = trim(adjustl(istring))//'/equiang/'
		else if (func.grid_type == 0) then
			if (func.rescale == 1) then
				istring = trim(adjustl(istring))//'/tan/'
			else if (func.rescale == 0) then
				istring = trim(adjustl(istring))//'/simple/'
			else if (func.rescale == 2) then
				istring = trim(adjustl(istring))//'/exp/'
			end if
		end if

		this.square = grid.square

		if(id == 0) then

! 			call this.histogram(16*grid.dim*grid.dim, 'datFiles/angle'//trim(istring)//'.dat', 'datFiles/angle_distribution'//trim(istring)//'.dat')
! 			call this.histogram(16*grid.dim*grid.dim, 'datFiles/dist'//trim(istring)//'.dat', 'datFiles/dist_distribution'//trim(istring)//'.dat')
! 			call this.histogram(4*grid.dim*grid.dim, 'datFiles/cell'//trim(istring)//'.dat', 'datFiles/cell_distribution'//trim(istring)//'.dat')

			open(9,file='datFiles/'//trim(istring)//'CFL.dat')
			open(11,file='datFiles/'//trim(istring)//'L1.dat')
			open(12,file='datFiles/'//trim(istring)//'L2.dat')
			open(13,file='datFiles/'//trim(istring)//'C.dat')
		end if

	end Subroutine



	Subroutine alloc(this, func)

		Class(diagnostic) :: this
		Class(f_var) :: func

		Allocate(this.CFL(func.first_x:func.last_x, func.first_y:func.last_y, 1:6))
		Allocate(this.square(1:2*func.dim, 1:2*func.dim))

	end Subroutine



	Subroutine deinit(this)
		Class(diagnostic) :: this
		Integer(4) ier, id

		if (Allocated(this.CFL)) Deallocate(this.CFL)
		call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)

		if(id == 0) then
			close(9)
			close(11)
			close(12)
			close(13)
		end if

	end Subroutine



	Subroutine Courant(this, metr, func, time)
		Class(diagnostic) :: this
		Class(f_var) :: func
		Class(metric) :: metr
		Integer(4), intent(in) :: time
		Integer(4) face, x, y, dim, ier, id
		Real(8) :: Courant_number, Courant_max, L1, L2, L1_prec, L2_prec, Linf, Linf_prec, square, F1, F1_prec, L1r, L2r, L1_precr, L2_precr, Linfr, Linf_precr, latlon1(1:2), r, pi, h, u0, alpha

		dim = this.dim;  Linf = 0d0;  Linf_prec = 0d0
		pi = 314159265358979323846d-20;  u0 = 2d0*pi/(12d0*24d0*60d0*60d0);  alpha = -pi/4d0

! 		if(time == 1) then
! 			this.center = (/0d0, 3d0*pi/2d0/)
! 		else
! 			latlon1(2) = u0*(dcos(this.center(1))*dcos(alpha) + dcos(this.center(2))*dsin(this.center(1))*dsin(alpha))*this.dt
! 			latlon1(1) = -u0*dsin(this.center(2))*dsin(alpha)*this.dt
! 		end if
! 		this.center(1) = latlon1(1) + this.center(1)
! 		this.center(2) = latlon1(2) + this.center(2)

! 		if(time == this.Tmax) this.center = (/0d0, 3d0*pi/2d0/)
! 		if(time == this.Tmax/2) this.center = (/0d0, pi/2d0/)
		! if(mod(time, 10) == 0 )print *, this.center

		!$OMP PARALLEL PRIVATE(y, x, face)
		!$OMP DO

		do face = 1, 6
			do y = func.ns_y, func.nf_y
				do x = func.ns_x, func.nf_x

					this.CFL(x, y, face) = dsqrt((func.u_con(x, y, face)*metr.G_sqr(x,y)*this.dt/(this.dh))**2 + (func.v_con(x, y, face)*metr.G_sqr(x,y)*this.dt/(this.dh))**2)

! 					r = dacos(dsin(this.center(1))*dsin(metr.latlon_c(1, x, y, face)) + dcos(this.center(1))*dcos(metr.latlon_c(1, x, y, face))*dcos(this.center(2) - metr.latlon_c(2, x, y, face)))
! 					func.starter(1, x, y, face) = (func.height/2d0)*(1d0 + dcos(3d0*pi*r))
! 					if (r >= 1d0/3d0) func.starter(1, x, y, face) = 0d0

					square = this.square(x, y)
					F1 = dsqrt((func.lon_vel(x, y, face) - func.starter(2, x, y, face))**2 + (func.lat_vel(x, y, face) - func.starter(3, x, y, face))**2)
					F1_prec = dsqrt((func.starter(2, x, y, face))**2 + (func.starter(3, x, y, face))**2)

! 					F1 = dsqrt((func.h_height(x, y, face) - func.starter(1, x, y, face))**2)
! 					F1_prec = dsqrt((func.starter(1, x, y, face))**2)

					! F1 = dsqrt((func.h_height(x, y, face) - h)**2)
					! F1_prec = dsqrt((h)**2)

					L1 = abs(F1)*square + L1
					L2 = F1*F1*square + L2

					L1_prec = abs(F1_prec)*square + L1_prec
					L2_prec = F1_prec*F1_prec*square + L2_prec

					if (Linf < F1) Linf = F1
					if (Linf_prec < F1_prec) Linf_prec = F1_prec

				end do
			end do
		end do

		!$OMP END DO
		!$OMP END PARALLEL

		Courant_number = MAXVAL(this.CFL)
		call MPI_Allreduce(Courant_number, Courant_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
		call MPI_Allreduce(L1, L1r, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
		call MPI_Allreduce(L2, L2r, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
		call MPI_Allreduce(Linf, Linfr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
		call MPI_Allreduce(L1_prec, L1_precr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
		call MPI_Allreduce(L2_prec, L2_precr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
		call MPI_Allreduce(Linf_prec, Linf_precr, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
		call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)

		if (id == 0) write(9, FMT = "(f40.6, f40.6)"),dble(time)*this.convert_time, Courant_max
		if (id == 0) write(11, FMT = "(f40.6, f40.6)"),dble(time)*this.convert_time, L1r/L1_precr
		if (id == 0) write(12, FMT = "(f40.6, f40.6)"),dble(time)*this.convert_time, dsqrt(L2r/L2_precr)
		if (id == 0) write(13, FMT = "(f40.6, f40.6)"),dble(time)*this.convert_time, Linfr/Linf_precr


	end Subroutine





	Subroutine hist_generation(this, N, filename_input, filename_output)

		Class(diagnostic) :: this
		Integer N
		character (len=*) filename_input, filename_output

		Integer, parameter :: hist_points=1000
		Real(8), allocatable :: dat(:)
		Real(8) max, min, max_val
		Integer(4) distribution(1:hist_points), i, index

		max=-1d20;min=1d20
		distribution = 0
		allocate(dat(1: N))
		open(20, file = filename_input)
		do i =1, N
			read(20,*) dat(i)
			if (dat(i)> max) then; max = dat(i); end if
			if (dat(i)< min) then; min = dat(i); end if
		end do
		close(20)

		do i = 1, N
			index = int(hist_points * (dat(i)-min) /(max-min) + 5d-1)
			distribution(index) = distribution(index) + 1
		end do

		max_val = maxval(distribution)

		open(21, file = filename_output)
		do i = 1, hist_points
			 if (distribution(i)>0) then
					write(21,*) (min + (max-min) * i /dble(hist_points)), distribution(i)
			 end if
		end do
		close(21)
		deallocate(dat)
	end Subroutine



end module