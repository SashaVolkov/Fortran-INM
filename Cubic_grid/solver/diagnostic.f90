module diagnostic_mod

	use grid_var, Only: g_var
	use metrics, Only: metric
	use func_var, Only: f_var
	use parallel_cubic, Only: parallel
	use mpi


	implicit none
	Private
	Public :: diagnostic


	Type diagnostic

		Real(8), Allocatable :: CFL(:,:,:)
		integer(4) Tmax, dim, step, flag
		real(8) convert_time, L10, L20, L_inf0


		CONTAINS
			Procedure, Public :: init => init
			Procedure, Private :: alloc => alloc
			Procedure, Public :: deinit => deinit
			Procedure, Public :: Courant => Courant
			Procedure, Public :: L_norm => L_norm
			Procedure, Public :: histogram => hist_generation

	End Type


CONTAINS



	subroutine init(this, grid, paral, Tmax, id)

		Class(diagnostic) :: this
		Class(g_var) :: grid
		Class(parallel) :: paral

		integer(4), intent(in) :: Tmax, id
		character(32) istring


		this.Tmax = Tmax;  this.dim = grid.dim;  this.step = grid.step
		this.convert_time = grid.dt/3600d0/24d0
		this.flag = 1

		write(istring, *) 2*this.dim

		call this.alloc(paral)

		if (grid.grid_type == 1) then
			istring = trim(adjustl(istring))//'/equiang/'
		else if (grid.grid_type == 0) then
			if (grid.rescale == 1) then
				istring = trim(adjustl(istring))//'/tan/'
			else if (grid.rescale == 0) then
				istring = trim(adjustl(istring))//'/simple/'
			else if (grid.rescale == 2) then
				istring = trim(adjustl(istring))//'/exp/'
			end if
		end if



		if(id == 0) then

! 			call this.histogram(16*grid.dim*grid.dim, 'datFiles/angle'//trim(istring)//'.dat', 'datFiles/angle_distribution'//trim(istring)//'.dat')
! 			call this.histogram(16*grid.dim*grid.dim, 'datFiles/dist'//trim(istring)//'.dat', 'datFiles/dist_distribution'//trim(istring)//'.dat')
! 			call this.histogram(4*grid.dim*grid.dim, 'datFiles/cell'//trim(istring)//'.dat', 'datFiles/cell_distribution'//trim(istring)//'.dat')

			open(9,file='datFiles/'//trim(istring)//'/CFL.dat')
			! open(11,file='datFiles/'//trim(istring)//'L1.dat')
			! open(12,file='datFiles/'//trim(istring)//'L2.dat')
			open(13,file='datFiles/'//trim(istring)//'/L_inf_cube.dat')

		end if

	end subroutine



	subroutine alloc(this, paral)

		Class(diagnostic) :: this
		Class(parallel) :: paral

		Allocate(this.CFL(paral.first_x:paral.last_x, paral.first_y:paral.last_y, 1:6))

	end subroutine



	subroutine deinit(this)
		Class(diagnostic) :: this
		integer(4) ier, id

		if (Allocated(this.CFL)) Deallocate(this.CFL)
		call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)

		if(id == 0) then
			close(9)
			! close(11)
			! close(12)
			close(13)
		end if

	end subroutine



	subroutine Courant(this, func, grid, metr, time)
		Class(diagnostic) :: this
		Class(f_var) :: func
		Class(g_var) :: grid
		Class(metric) :: metr
		integer(4), intent(in) :: time
		integer(4) face, x, y, dim, ier, id
		real(8) :: Courant_number, Courant_max

		dim = this.dim

		do face = 1, 6
			do y = func.ns_y, func.nf_y
				do x = func.ns_x, func.nf_x

					this.CFL(x, y, face) = abs(func.u_cov(x, y, face)*grid.dt/(grid.delta_on_cube)) +&
					 abs(func.v_cov(x, y, face)*grid.dt/(grid.delta_on_cube))

				end do
			end do
		end do

		Courant_number = MAXVAL(this.CFL)
		call MPI_Allreduce(Courant_number, Courant_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
		call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)

		if (id == 0) write(9, FMT = "(f40.6, f40.6)"),dble(time)*this.convert_time, Courant_max

	end subroutine


	subroutine L_norm(this, func, grid, time)
		Class(diagnostic) :: this
		Class(g_var) :: grid
		Class(f_var) :: func
		integer(4), intent(in) :: time

		integer(4) face, x, y, id, ier
		real(8) L1, L2, L1_all, L2_all, L_inf, L_inf_all, F1, F2, square

		L_inf = MAXVAL(abs(func.h_height))
		call MPI_Allreduce(L_inf, L_inf_all, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)

		if (id == 0 ) then
			if ( this.flag == 1 ) then
				write(13, FMT = *), 100
				this.flag = 0
			end if
			write(13, FMT = *), abs(L_inf_all)
		end if


	end subroutine



	subroutine hist_generation(this, N, filename_input, filename_output)

		Class(diagnostic) :: this
		integer N
		character (len=*) filename_input, filename_output

		integer, parameter :: hist_points=1000
		real(8), allocatable :: dat(:)
		real(8) max, min, max_val
		integer(4) distribution(1:hist_points), i, index

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
	end subroutine



end module