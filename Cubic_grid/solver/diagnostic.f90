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
		Real(8) radius, pi


		CONTAINS
			Procedure, Public :: init => init
			Procedure, Private :: alloc => alloc
			Procedure, Public :: deinit => deinit
			Procedure, Public :: CFL => CFL
			Procedure, Private :: spherical_triangle => spherical_triangle
			Procedure, Private :: spherical_triangle_area => spherical_triangle_area

	End Type


CONTAINS



	subroutine init(this, grid, Tmax, rescale, pi)

		Class(diagnostic) :: this
		Class(g_var) :: grid
		integer(4), intent(in) :: Tmax, rescale
		real(8), intent(in) :: pi

		this.Tmax = Tmax;  this.dim = grid.dim;  this.step = grid.step
		this.radius = grid.r_sphere;  this.pi = pi

		call this.alloc()
		if (rescale == 1) then
			open(9,file='/home/sasha/Fortran/Cubic_grid/solver/datFiles/CFL_x.dat')
			open(10,file='/home/sasha/Fortran/Cubic_grid/solver/datFiles/CFL_y.dat')
		else
			open(9,file='/home/sasha/Fortran/Cubic_grid/solver/datFiles/CFL_x_simple.dat')
			open(10,file='/home/sasha/Fortran/Cubic_grid/solver/datFiles/CFL_y_simple.dat')
		end if

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


	!!!!Wiki article "Решение треугольников"
	subroutine spherical_triangle(this, a, b, c, alpha, beta, gamma)

		Class(diagnostic) :: this
		real(8), intent(in) :: a, b, c  ! angles between radiuses
		real(8), intent(out) :: alpha, beta, gamma  ! angles of spherical triangle

		alpha = dacos( ( dcos(a) - dcos(b)*dcos(c) )/( dsin(b)*dsin(c) ) )
		beta = dacos( ( dcos(b) - dcos(a)*dcos(c) )/( dsin(a)*dsin(c) ) )
		gamma = dacos( ( dcos(c) - dcos(b)*dcos(a) )/( dsin(b)*dsin(a) ) )

	end subroutine



		subroutine spherical_triangle_area(this, alpha, beta, gamma, area)

		Class(diagnostic) :: this
		real(8), intent(in) :: alpha, beta, gamma  ! angles of spherical triangle
		real(8), intent(out) :: area
		real(8) eps

		eps = alpha + beta + gamma - this.pi
		area = this.radius * this.radius * eps

	end subroutine



end module