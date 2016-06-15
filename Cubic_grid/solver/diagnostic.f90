module diagnostic_mod

	use grid_var, Only: g_var
	use func_var, Only: f_var

	implicit none
	Private
	Public :: diagnostic

	Type diagnostic

		Real(8), Allocatable :: CFL(:,:,:)
		integer(4) Tmax, dim, step
		real(8) convert_time


		CONTAINS
			Procedure, Public :: init => init
			Procedure, Private :: alloc => alloc
			Procedure, Public :: deinit => deinit
			Procedure, Public :: Courant => Courant
			Procedure, Public :: L_norm => L_norm
			Procedure, Public :: histogram => hist_generation

	End Type


CONTAINS



	subroutine init(this, grid, Tmax, rescale)

		Class(diagnostic) :: this
		Class(g_var) :: grid
		integer(4), intent(in) :: Tmax, rescale
		character(8) istring


		this.Tmax = Tmax;  this.dim = grid.dim;  this.step = grid.step
		this.convert_time = grid.dt/3600d0/24d0

		call this.alloc()

		if (rescale == 1) then
			istring = '_tan'
		else if (rescale == 0) then
			istring = '_simple'
		else if (rescale == 2) then
			istring = '_exp'
		end if

		call this.histogram(16*grid.dim*grid.dim, 'datFiles/angle'//trim(istring)//'.dat', 'datFiles/angle_distribution'//trim(istring)//'.dat')
		call this.histogram(16*grid.dim*grid.dim, 'datFiles/dist'//trim(istring)//'.dat', 'datFiles/dist_distribution'//trim(istring)//'.dat')
		call this.histogram(4*grid.dim*grid.dim, 'datFiles/cell'//trim(istring)//'.dat', 'datFiles/cell_distribution'//trim(istring)//'.dat')

		open(9,file='datFiles/CFL'//trim(istring)//'.dat')
		open(11,file='datFiles/L1'//trim(istring)//'.dat')
		open(12,file='datFiles/L2'//trim(istring)//'.dat')
		open(13,file='datFiles/L_inf'//trim(istring)//'.dat')

	end subroutine



	subroutine alloc(this)

		Class(diagnostic) :: this

		Allocate(this.CFL(-this.dim:this.dim, -this.dim:this.dim, 1:6))

	end subroutine



	subroutine deinit(this)
		Class(diagnostic) :: this

		if (Allocated(this.CFL)) Deallocate(this.CFL)
		close(9)

	end subroutine



	subroutine Courant(this, func, grid, time)
		Class(diagnostic) :: this
		Class(f_var) :: func(1:6)
		Class(g_var) :: grid
		integer(4), intent(in) :: time
		integer(4) face, x, y, dim

		dim = this.dim

		do face = 1, 6
			do y = -dim, dim
				do x = -dim, dim

					this.CFL(x, y, face) = abs(func(face).u_vel(x, y)*grid.dt/grid.h_dist(2, 1, y, x)) +&
					 abs(func(face).v_vel(x, y)*grid.dt/grid.h_dist(3, 1, y, x))

				end do
			end do
		end do

			write(9, FMT = "(f14.6, f10.4)"),dble(time)*this.convert_time, MAXVAL(this.CFL)


! 		print '(" CFL x = ", f6.4, "   CFL y = ", f6.4)', MAXVAL(this.CFL_x), MAXVAL(this.CFL_y)

	end subroutine



	subroutine L_norm(this, func, grid, time)
		Class(diagnostic) :: this
		Class(g_var) :: grid
		Class(f_var) :: func(1:6)
! 		real(8), intent(in) :: func(-this.dim:this.dim, -this.dim:this.dim, 1:6)
		integer(4), intent(in) :: time
		integer(4) face, x, y, dim
		real(8) L1, L2, F1, F2, L_inf

		L1 = 0;  dim = this.dim;  L_inf = 0

		do face = 1, 6
			do y = -dim, dim
				do x = -dim, dim

					F1 = func(face).h_height(x, y) + func(face).h_height(x+1, y) + func(face).h_height(x, y+1)
					F2 = func(face).h_height(x+1, y+1) + func(face).h_height(x+1, y) + func(face).h_height(x, y+1)
					L1 = abs(F1)*grid.triangle_area(1, x, y) + abs(F2)*grid.triangle_area(2, x, y)
					L2 = F1*F1*grid.triangle_area(1, x, y) + F2*F2*grid.triangle_area(2, x, y)

					if ( L_inf < MAXVAL(abs(func(face).h_height)) ) then
						L_inf = MAXVAL(abs(func(face).h_height))
					end if


				end do
			end do
		end do

		L2 = dsqrt(L2)

		write(11, FMT = "(f14.6, f10.4)"),time*this.convert_time, L1
		write(12, FMT = "(f14.6, f10.4)"),time*this.convert_time, L2
		write(13, FMT = "(f14.6, f10.4)"),time*this.convert_time, L_inf

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