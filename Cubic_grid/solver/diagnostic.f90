module diagnostic_mod

	use grid_var, Only: g_var
	use func_var, Only: f_var

	implicit none
	Private
	Public :: diagnostic

	Type diagnostic

		Real(8), Allocatable :: CFL_x(:,:,:)
		Real(8), Allocatable :: CFL_y(:,:,:)
		integer(4) Tmax, dim, step


		CONTAINS
			Procedure, Public :: init => init
			Procedure, Private :: alloc => alloc
			Procedure, Public :: deinit => deinit
			Procedure, Public :: CFL => CFL
			Procedure, Public :: L_norm1 => L_norm1

	End Type


CONTAINS



	subroutine init(this, grid, Tmax, rescale)

		Class(diagnostic) :: this
		Class(g_var) :: grid
		integer(4), intent(in) :: Tmax, rescale
		character(8) istring


		this.Tmax = Tmax;  this.dim = grid.dim;  this.step = grid.step

		call this.alloc()

		if (rescale == 1) then
			istring = '_tan'
		else if (rescale == 0) then
			istring = '_simple'
		else if (rescale == 2) then
			istring = '_exp'
		end if

		open(9,file='/home/sasha/Fortran/Cubic_grid/solver/datFiles/CFL_x'//trim(istring)//'.dat')
		open(10,file='/home/sasha/Fortran/Cubic_grid/solver/datFiles/CFL_y'//trim(istring)//'.dat')
		open(11,file='/home/sasha/Fortran/Cubic_grid/solver/datFiles/L1'//trim(istring)//'.dat')
		open(12,file='/home/sasha/Fortran/Cubic_grid/solver/datFiles/L2'//trim(istring)//'.dat')
		open(13,file='/home/sasha/Fortran/Cubic_grid/solver/datFiles/L_inf'//trim(istring)//'.dat')

	end subroutine



	subroutine alloc(this)

		Class(diagnostic) :: this

		Allocate(this.CFL_x(-this.dim:this.dim, -this.dim:this.dim, 1:6))
		Allocate(this.CFL_y(-this.dim:this.dim, -this.dim:this.dim, 1:6))

	end subroutine



	subroutine deinit(this)
		Class(diagnostic) :: this

		if (Allocated(this.CFL_x)) Deallocate(this.CFL_x)
		if (Allocated(this.CFL_y)) Deallocate(this.CFL_y)
		close(9)
		close(10)

	end subroutine



	subroutine CFL(this, func, grid, time)
		Class(diagnostic) :: this
		Class(f_var) :: func
		Class(g_var) :: grid
		integer(4), intent(in) :: time
		integer(4) face_idx, x, y, dim

		dim = this.dim

		do face_idx = 1, 6
			do y = -dim, dim
				do x = -dim, dim

					this.CFL_x(x, y, face_idx) = abs(func.u_vel(x, y, face_idx)*grid.dt/grid.h_dist(2, y, x))
					this.CFL_y(x, y, face_idx) = abs(func.v_vel(x, y, face_idx)*grid.dt/grid.h_dist(3, y, x))

				end do
			end do
		end do

			write(9, FMT = "(I14, f10.4)"),time, MAXVAL(this.CFL_x)
			write(10, FMT = "(I14, f10.4)"),time, MAXVAL(this.CFL_y)


! 		print '(" CFL x = ", f6.4, "   CFL y = ", f6.4)', MAXVAL(this.CFL_x), MAXVAL(this.CFL_y)

	end subroutine



	subroutine L_norm1(this, func, grid, time)
		Class(diagnostic) :: this
		Class(g_var) :: grid
		real(8), intent(in) :: func(-this.dim:this.dim, -this.dim:this.dim, 1:6)
		integer(4), intent(in) :: time
		integer(4) face_idx, x, y, dim
		real(8) L1, L2, F1, F2

		L1 = 0;  dim = this.dim

		do face_idx = 1, 6
			do y = -dim, dim
				do x = -dim, dim

					F1 = func(x, y, face_idx) + func(x+1, y, face_idx) + func(x, y+1, face_idx)
					F2 = func(x+1, y+1, face_idx) + func(x+1, y, face_idx) + func(x, y+1, face_idx)
					L1 = abs(F1)*grid.triangle_area(1, x, y) + abs(F2)*grid.triangle_area(2, x, y)
					L2 = F1*F1*grid.triangle_area(1, x, y) + F2*F2*grid.triangle_area(2, x, y)

				end do
			end do
		end do

		L2 = dsqrt(L2)

		write(11, FMT = "(I14, f10.4)"),time, L1
		write(12, FMT = "(I14, f10.4)"),time, L2
		write(13, FMT = "(I14, f10.4)"),time, MAXVAL(abs(func))

	end subroutine


end module