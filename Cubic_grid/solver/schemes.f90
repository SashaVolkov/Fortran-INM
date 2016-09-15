module schemes

	use grid_var, Only: g_var
	use derivatives, Only: der
	use metrics, Only: metric
	use func_var, Only: f_var
	use interpolation, Only: interp
	use parallel_cubic, Only: parallel
	use messenger, Only: message
	use mpi

	implicit none


	Private
	Public :: schema

	Type schema

		integer(4) first_x, first_y, last_x, last_y

		Real(8), Allocatable :: ku_cov(:, :, :, :)
		Real(8), Allocatable :: kv_cov(:, :, :, :)
		Real(8), Allocatable :: ku_con(:, :, :, :)
		Real(8), Allocatable :: kv_con(:, :, :, :)
		Real(8), Allocatable :: kh(:, :, :, :)

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

	this.first_x = f_x;  this.first_y = f_y
	this.last_x = l_x;  this.last_y = l_y

	Allocate(this.ku_cov(f_x: l_x, f_y : l_y, 6, 0:4))
	Allocate(this.kv_cov(f_x: l_x, f_y : l_y, 6, 0:4))
	Allocate(this.ku_con(f_x: l_x, f_y : l_y, 6, 0:4))
	Allocate(this.kv_con(f_x: l_x, f_y : l_y, 6, 0:4))
	Allocate(this.kh(f_x: l_x, f_y : l_y, 6, 0:4))


End Subroutine



subroutine Linear(this, var, var_pr, grid, metr)

	Class(schema) :: this
	Class(f_var) :: var, var_pr
	Class(g_var) :: grid
	Class(metric) :: metr
	Type(der) :: d

	real(8) g, height, dt, partial, temp1(-2:2), temp2(-2:2), div, h
	integer(4) face, x, y, dim, order

	g = grid.g;  height = var_pr.height;  dim = var_pr.dim
	dt = grid.dt;  order = 2


	do face = 1, 6
		do y = var.ns_y, var.nf_y
			do x = var.ns_x, var.nf_x

				h = grid.delta_on_cube
				temp1(:) = var_pr.h_height(x-order:x+order, y, face)
				partial = d.partial_c4(temp1, h)
				var.u_cov(x, y, face) = var_pr.u_cov(x, y, face) - dt*g*partial

				temp1(:) = var_pr.h_height(x, y-order:y+order, face)
				partial = d.partial_c4(temp1, h)
				var.v_cov(x, y, face) = var_pr.v_cov(x, y, face) - dt*g*partial

				temp1(:) = var_pr.u_con(x-order:x+order, y, face)
				temp2(:) = var_pr.v_con(x, y-order:y+order, face)
				div = d.div(metr, temp1, temp2, h, x, y, order)
				var.h_height(x, y, face) = var_pr.h_height(x, y, face) - dt*height*div
			end do
		end do
	end do


end subroutine




Subroutine RungeKutta(this, var, var_pr, grid, metr, inter, paral, msg)

	Class(schema) :: this
	Class(f_var) :: var, var_pr
	Class(g_var) :: grid
	Class(metric) :: metr
	Class(interp) :: inter
	Class(parallel) :: paral
	Class(message) :: msg

	integer(4) face, x, y, dim, i, j, stat, ns_x, ns_y, nf_x, nf_y, ier, iteration

	dim = var_pr.dim
	ns_x = var.ns_x;  ns_y = var.ns_y
	nf_x = var.nf_x;  nf_y = var.nf_y

	call var_pr.cov_to_con(metr)


	this.ku_cov(:, :, :, 0) = var_pr.u_cov(:, :, :)
	this.kv_cov(:, :, :, 0) = var_pr.v_cov(:, :, :)
	this.ku_con(:, :, :, 0) = var_pr.u_con(:, :, :)
	this.kv_con(:, :, :, 0) = var_pr.v_con(:, :, :)
	this.kh(:, :, :, 0) = var_pr.h_height(:, :, :)

		do iteration = 1, 4
			call this.FRunge(grid, metr, var_pr, iteration)
			var_pr.u_cov(:, :, :) = this.ku_cov(:, :, :, iteration)
			var_pr.v_cov(:, :, :) = this.kv_cov(:, :, :, iteration)
			var_pr.h_height(:, :, :) = this.kh(:, :, :, iteration)
			call var_pr.equal(var_pr, metr)
			call msg.msg(var_pr, paral)
			call var_pr.interpolate(inter, metr)
			this.ku_cov(:, :, :, iteration) = var_pr.u_cov(:, :, :)
			this.kv_cov(:, :, :, iteration) = var_pr.v_cov(:, :, :)
			this.ku_con(:, :, :, iteration) = var_pr.u_con(:, :, :)
			this.kv_con(:, :, :, iteration) = var_pr.v_con(:, :, :)
			this.kh(:, :, :, iteration) = var_pr.h_height(:, :, :)
		end do


	do face = 1, 6
		do y = ns_y, nf_y
			do x= ns_x, nf_x
var.u_cov(x, y, face) = this.ku_cov(x, y, face, 0) + (this.ku_cov(x, y, face, 1) + 2.0*this.ku_cov(x, y, face, 2) + 2.0*this.ku_cov(x, y, face, 3) + this.ku_cov(x, y, face, 4))/6.0
var.v_cov(x, y, face) = this.kv_cov(x, y, face, 0) + (this.kv_cov(x, y, face, 1) + 2.0*this.kv_cov(x, y, face, 2) + 2.0*this.kv_cov(x, y, face, 3) + this.kv_cov(x, y, face, 4))/6.0
var.h_height(x, y, face) = this.kh(x, y, face, 0) + (this.kh(x, y, face, 1) + 2.0*this.kh(x, y, face, 2) + 2.0*this.kh(x, y, face, 3) + this.kh(x, y, face, 4))/6.0
			end do
		end do
	end do

End Subroutine


Subroutine FRunge(this, grid, metr, var, i)
	Class(schema) :: this
	Class(f_var) :: var
	Class(g_var) :: grid
	Class(metric) :: metr
	Type(der) :: d

	integer(4), intent(in) :: i
	real(8) g, height, dt, partial, temp1(-2:2), temp2(-2:2), coef(0:3), div, h
	integer(4) x,y, face

	coef(0) = 0d0;  coef(1) = 0d5;  coef(2) = 0d5;  coef(3) = 1d0;

	dt = grid.dt;  g = grid.g; height = var.height
	h = grid.delta_on_cube

	do face = 1, 6
		do y = var.ns_y, var.nf_y
			do x = var.ns_x, var.nf_x

				temp1(:) = this.kh(x-2:x+2, y, face, 0) + coef(i-1)*this.kh(x-2:x+2, y, face, i-1)
				partial = d.partial_c4(temp1, h)
				this.ku_cov(x, y, face, i) = - dt*g*partial

				temp2(:) = this.kh(x, y-2:y+2, face, 0) + coef(i-1)*this.kh(x, y-2:y+2, face, i-1)
				partial = d.partial_c4(temp2, h)
				this.kv_cov(x, y, face, i) =  - dt*g*partial

				temp1(:) = this.ku_con(x-2:x+2, y, face, 0) + coef(i-1)*this.ku_con(x-2:x+2, y, face, i-1)
				temp2(:) = this.kv_con(x, y-2:y+2, face, 0) + coef(i-1)*this.kv_con(x, y-2:y+2, face, i-1)
				div = d.div(metr, temp1(:), temp2(:), h, x, y, 2)
				this.kh(x, y, face, i) = - dt*height*div

			end do
		end do
	end do

end Subroutine




Subroutine deinit(this)
	Class(schema) :: this

	if (Allocated(this.ku_cov)) Deallocate(this.ku_cov)
	if (Allocated(this.kv_cov)) Deallocate(this.kv_cov)
	if (Allocated(this.ku_con)) Deallocate(this.ku_con)
	if (Allocated(this.kv_con)) Deallocate(this.kv_con)
	if (Allocated(this.kh)) Deallocate(this.kh)

End Subroutine



end module
