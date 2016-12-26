module diagnostic_mod

	use func_var, Only: f_var
	use mpi
	use omp_lib


	implicit none
	Private
	Public :: diagnostic


	Type diagnostic

		Real(8), Allocatable :: CFL(:,:,:)
		Integer(4) Tmax, dim, step, flag
		Real(8) convert_time, L10, L20, L_inf0, dh, dt


		CONTAINS
			Procedure, Public :: init => init
			Procedure, Private :: alloc => alloc
			Procedure, Public :: deinit => deinit
			Procedure, Public :: Courant => Courant
			Procedure, Public :: histogram => hist_generation

	End Type


CONTAINS



	Subroutine init(this, func, Tmax, id)

		Class(diagnostic) :: this
		Class(f_var) :: func

		Integer(4), intent(in) :: Tmax, id
		character(32) istring, istring1

		this.Tmax = Tmax;  this.dim = func.dim;  this.step = func.step-1
		this.dt = func.dt;  this.dh = func.delta_on_cube
		this.convert_time = this.dt/(644d0*4000d0)
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



		if(id == 0) then

! 			call this.histogram(16*grid.dim*grid.dim, 'datFiles/angle'//trim(istring)//'.dat', 'datFiles/angle_distribution'//trim(istring)//'.dat')
! 			call this.histogram(16*grid.dim*grid.dim, 'datFiles/dist'//trim(istring)//'.dat', 'datFiles/dist_distribution'//trim(istring)//'.dat')
! 			call this.histogram(4*grid.dim*grid.dim, 'datFiles/cell'//trim(istring)//'.dat', 'datFiles/cell_distribution'//trim(istring)//'.dat')

			open(9,file='datFiles/'//trim(istring)//'CFL.dat')

		end if

	end Subroutine



	Subroutine alloc(this, func)

		Class(diagnostic) :: this
		Class(f_var) :: func

		Allocate(this.CFL(func.first_x:func.last_x, func.first_y:func.last_y, 1:6))

	end Subroutine



	Subroutine deinit(this)
		Class(diagnostic) :: this
		Integer(4) ier, id

		if (Allocated(this.CFL)) Deallocate(this.CFL)
		call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)

		if(id == 0) then
			close(9)
			! close(11)
			! close(12)
			! close(13)
		end if

	end Subroutine



	Subroutine Courant(this, func, time)
		Class(diagnostic) :: this
		Class(f_var) :: func
		Integer(4), intent(in) :: time
		Integer(4) face, x, y, dim, ier, id
		Real(8) :: Courant_number, Courant_max

		dim = this.dim

		!$OMP PARALLEL PRIVATE(y, x, face)
		!$OMP DO

		do face = 1, 6
			do y = func.ns_y, func.nf_y
				do x = func.ns_x, func.nf_x

					this.CFL(x, y, face) = abs(func.u_cov(x, y, face)*this.dt/(this.dh)) +&
					 abs(func.v_cov(x, y, face)*this.dt/(this.dh))

				end do
			end do
		end do

		!$OMP END DO
		!$OMP END PARALLEL

		Courant_number = MAXVAL(this.CFL)
		call MPI_Allreduce(Courant_number, Courant_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
		call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)

		if (id == 0) write(9, FMT = "(f40.6, f40.6)"),dble(time)*this.convert_time, Courant_max

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