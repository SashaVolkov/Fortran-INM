module schemes

	use grid_var, Only: g_var
	use func_var, Only: f_var
	use mpi

	implicit none


	Private
	Public :: schema

	Type schema

		Real(8), Allocatable :: ku(:, :, :)
		Real(8), Allocatable :: kv(:, :, :)
		Real(8), Allocatable :: kh(:, :, :)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: Linear => Linear
		Procedure, Public ::  RungeKutta=> RungeKutta
		Procedure, Private ::  FRunge=> FRunge
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS


Subroutine init(this, f, g)

	Class(g_var) :: g
	Class(f_var) :: f
	Class(schema) :: this

	integer(4) f_x, f_y, l_x, l_y

	f_x = f.first_x;  f_y = f.first_y
	l_x = f.last_x;  l_y = f.last_y

	Allocate(this.ku(f_x: l_x, f_y : l_y, 0 : 4))
	Allocate(this.kv(f_x: l_x, f_y : l_y, 0 : 4))
	Allocate(this.kh(f_x: l_x, f_y : l_y, 0 : 4))


End Subroutine



subroutine Linear(this, var, var_pr, grid)

	Class(schema) :: this
	Class(f_var) :: var, var_pr
	Class(g_var) :: grid

	real(8) g, height, dt, partial(1:2), temp1(-1:1), temp2(-1:1), div
	integer(4) face, x, y, dim, i, j, stat, ier

	g = grid.g;  height = var_pr.height;  dim = var_pr.dim
	dt = grid.dt

	do face = 2, 2
	do y = var.ns_y, var.nf_y
		do x = var.ns_x, var.nf_x

				partial(1) = grid.partial_c2_x(var_pr.h_height(x-1:x+1, y, face), x, y)
				var.x_vel(x, y, face) = var_pr.x_vel(x, y, face) - dt*g*partial(1)

				temp1(:) = var_pr.h_height(x, y-1:y+1, face)
				partial(1) = grid.partial_c2_y(temp1, x, y)
				var.y_vel(x, y, face) = var_pr.y_vel(x, y, face) - dt*g*partial(1)

				temp1(:) = var_pr.x_vel(x-1:x+1, y, face)
				temp2(:) = var_pr.y_vel(x, y-1:y+1, face)
				div = grid.div_2(temp1, temp2, x, y)
				var.h_height(x, y, face) = var_pr.h_height(x, y, face) - dt*height*div


			end do
		end do
	end do

		call var_pr.equal(var)

end subroutine




Subroutine RungeKutta(this, var, var_pr, grid)

	Class(schema) :: this
	Class(f_var) :: var, var_pr
	Class(g_var) :: grid

	integer(4) face, x, y, dim, i, j, stat, ns_x, ns_y, nf_x, nf_y, ier, iteration

	dim = var_pr.dim
	ns_x = var.ns_x;  ns_y = var.ns_y
	nf_x = var.nf_x;  nf_y = var.nf_y


	do face = 2, 2

	this.ku(:, :, 0) = var_pr.x_vel(:, :, face)
	this.kv(:, :, 0) = var_pr.y_vel(:, :, face)
	this.kh(:, :, 0) = var_pr.h_height(:, :, face)

		do iteration = 1, 4
			call this.FRunge(grid, var_pr, face, iteration)
		end do


		do y = ns_y, nf_y
			do x= ns_x, nf_x
var.x_vel(x, y, face) = var_pr.x_vel(x, y, face) + (this.ku(x, y, 1) + 2.0*this.ku(x, y, 2) + 2.0*this.ku(x, y, 3) + this.ku(x, y, 4))/6.0
var.y_vel(x, y, face) = var_pr.y_vel(x, y, face) + (this.kv(x, y, 1) + 2.0*this.kv(x, y, 2) + 2.0*this.kv(x, y, 3) + this.kv(x, y, 4))/6.0
var.h_height(x, y, face) = var_pr.h_height(x, y, face) + (this.kh(x, y, 1) + 2.0*this.kh(x, y, 2) + 2.0*this.kh(x, y, 3) + this.kh(x, y, 4))/6.0
			end do
		end do
	end do

	call var_pr.equal(var)

End Subroutine


Subroutine FRunge(this, grid, var, face, i)
	Class(schema) :: this
	Class(f_var) :: var
	Class(g_var) :: grid

	integer(4), intent(in) :: face, i
	real(8) g, height, dt, partial(1:2), temp1(-2:2), temp2(-2:2), coef(0:3), div
	integer(4) x,y

	coef(0) = 0d0;  coef(1) = 0d5;  coef(2) = 0d5;  coef(3) = 1d0;

	dt = grid.dt;  g = grid.g; height = var.height

	do y = var.ns_y, var.nf_y
		do x = var.ns_x, var.nf_x

			temp1(:) = this.kh(x-2:x+2, y, 0) + coef(i-1)*this.kh(x-2:x+2, y, i-1)
			partial(1) = grid.partial_c4_x(temp1, x, y)
			this.ku(x, y, i) = - dt*g*partial(1)

			temp2(:) = this.kh(x, y-2:y+2, 0) + coef(i-1)*this.kh(x, y-2:y+2, i-1)
			partial(1) = grid.partial_c4_y(temp2, x, y)
			this.kv(x, y, i) =  - dt*g*partial(1)

			temp1(:) = this.ku(x-2:x+2, y, 0) + coef(i-1)*this.ku(x-2:x+2, y, i-1)
			temp2(:) = this.kv(x, y-2:y+2, 0) + coef(i-1)*this.kv(x, y-2:y+2, i-1)
			div = grid.div_4(temp1(:), temp2(:), x, y)
			this.kh(x, y, i) = - dt*height*div

		end do
	end do

end Subroutine




Subroutine deinit(this)
	Class(schema) :: this

	if (Allocated(this.ku)) Deallocate(this.ku)
	if (Allocated(this.kv)) Deallocate(this.kv)
	if (Allocated(this.kh)) Deallocate(this.kh)

End Subroutine




end module
