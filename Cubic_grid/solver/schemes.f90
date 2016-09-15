module schemes

	use grid_var, Only: g_var
	use derivatives, Only: der
	use metrics, Only: metric
	use func_var, Only: f_var
	use mpi

	implicit none


	Private
	Public :: schema

	Type schema

		integer(4) first_x, first_y, last_x, last_y

		Real(8), Allocatable :: ku(:, :, :)
		Real(8), Allocatable :: kv(:, :, :)
		Real(8), Allocatable :: ku_con(:, :, :)
		Real(8), Allocatable :: kv_con(:, :, :)
		Real(8), Allocatable :: kh(:, :, :)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: Linear => Linear
		Procedure, Public ::  RungeKutta=> RungeKutta
		Procedure, Private ::  FRunge=> FRunge
		Procedure, Public :: deinit => deinit
		Procedure, Public :: cov_to_con => cov_to_con
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

	Allocate(this.ku(f_x: l_x, f_y : l_y, 0:4))
	Allocate(this.kv(f_x: l_x, f_y : l_y, 0:4))
	Allocate(this.ku_con(f_x: l_x, f_y : l_y, 0:4))
	Allocate(this.kv_con(f_x: l_x, f_y : l_y, 0:4))
	Allocate(this.kh(f_x: l_x, f_y : l_y, 0:4))


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




Subroutine RungeKutta(this, var, var_pr, grid, metr)

	Class(schema) :: this
	Class(f_var) :: var, var_pr
	Class(g_var) :: grid
	Class(metric) :: metr

	integer(4) face, x, y, dim, i, j, stat, ns_x, ns_y, nf_x, nf_y, ier, iteration

	dim = var_pr.dim
	ns_x = var.ns_x;  ns_y = var.ns_y
	nf_x = var.nf_x;  nf_y = var.nf_y

	call var_pr.cov_to_con(metr)

	do face = 1, 6

	this.ku(:, :, 0) = var_pr.u_cov(:, :, face)
	this.kv(:, :, 0) = var_pr.v_cov(:, :, face)
	this.ku_con(:, :, 0) = var_pr.u_con(:, :, face)
	this.kv_con(:, :, 0) = var_pr.v_con(:, :, face)
	this.kh(:, :, 0) = var_pr.h_height(:, :, face)

		do iteration = 1, 4
			call this.FRunge(grid, metr, var_pr, face, iteration)
			call this.cov_to_con(metr, i)
		end do


		do y = ns_y, nf_y
			do x= ns_x, nf_x
var.u_cov(x, y, face) = var_pr.u_cov(x, y, face) + (this.ku(x, y, 1) + 2.0*this.ku(x, y, 2) + 2.0*this.ku(x, y, 3) + this.ku(x, y, 4))/6.0
var.v_cov(x, y, face) = var_pr.v_cov(x, y, face) + (this.kv(x, y, 1) + 2.0*this.kv(x, y, 2) + 2.0*this.kv(x, y, 3) + this.kv(x, y, 4))/6.0
var.h_height(x, y, face) = var_pr.h_height(x, y, face) + (this.kh(x, y, 1) + 2.0*this.kh(x, y, 2) + 2.0*this.kh(x, y, 3) + this.kh(x, y, 4))/6.0
			end do
		end do
	end do

End Subroutine


Subroutine FRunge(this, grid, metr, var, face, i)
	Class(schema) :: this
	Class(f_var) :: var
	Class(g_var) :: grid
	Class(metric) :: metr
	Type(der) :: d

	integer(4), intent(in) :: face, i
	real(8) g, height, dt, partial, temp1(-2:2), temp2(-2:2), coef(0:3), div, h
	integer(4) x,y

	coef(0) = 0d0;  coef(1) = 0d5;  coef(2) = 0d5;  coef(3) = 1d0;

	dt = grid.dt;  g = grid.g; height = var.height
	h = grid.delta_on_cube

	do y = var.ns_y, var.nf_y
		do x = var.ns_x, var.nf_x

			temp1(:) = this.kh(x-2:x+2, y, 0) + coef(i-1)*this.kh(x-2:x+2, y, i-1)
			partial = d.partial_c4(temp1, h)
			this.ku(x, y, i) = - dt*g*partial

			temp2(:) = this.kh(x, y-2:y+2, 0) + coef(i-1)*this.kh(x, y-2:y+2, i-1)
			partial = d.partial_c4(temp2, h)
			this.kv(x, y, i) =  - dt*g*partial

			temp1(:) = this.ku_con(x-2:x+2, y, 0) + coef(i-1)*this.ku_con(x-2:x+2, y, i-1)
			temp2(:) = this.kv_con(x, y-2:y+2, 0) + coef(i-1)*this.kv_con(x, y-2:y+2, i-1)
			div = d.div(metr, temp1(:), temp2(:), h, x, y, 2)
			this.kh(x, y, i) = - dt*height*div



		end do
	end do

end Subroutine




Subroutine deinit(this)
	Class(schema) :: this

	if (Allocated(this.ku)) Deallocate(this.ku)
	if (Allocated(this.kv)) Deallocate(this.kv)
	if (Allocated(this.ku_con)) Deallocate(this.ku_con)
	if (Allocated(this.kv_con)) Deallocate(this.kv_con)
	if (Allocated(this.kh)) Deallocate(this.kh)

End Subroutine



	subroutine cov_to_con(this, metr, i)
		Class(schema) :: this
		Class(metric) :: metr
		integer(4), intent(in) :: i
		Real(8) :: vel_x_contr, vel_y_contr
		Integer(4) :: x, y, face

			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x

this.ku_con(x, y, i) = metr.G_inverse(1, 1, x, y) * this.ku(x, y, i) + metr.G_inverse(1, 2, x, y) * this.kv(x, y, i)
this.kv_con(x, y, i) = metr.G_inverse(2, 2, x, y) * this.kv(x, y, i) + metr.G_inverse(2, 1, x, y) * this.ku(x, y, i)

				end do
			end do

	end subroutine




end module
