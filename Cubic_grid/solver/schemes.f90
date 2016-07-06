module schemes

	use grid_var, Only: g_var
	use func_var, Only: f_var
	use mpi

	implicit none


	Private
	Public :: schema

	Type schema

		Real(8), Allocatable :: ku1(:,:)
		Real(8), Allocatable :: ku2(:,:)
		Real(8), Allocatable :: ku3(:,:)
		Real(8), Allocatable :: ku4(:,:)

		Real(8), Allocatable :: kv1(:,:)
		Real(8), Allocatable :: kv2(:,:)
		Real(8), Allocatable :: kv3(:,:)
		Real(8), Allocatable :: kv4(:,:)

		Real(8), Allocatable :: kh1(:,:)
		Real(8), Allocatable :: kh2(:,:)
		Real(8), Allocatable :: kh3(:,:)
		Real(8), Allocatable :: kh4(:,:)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: Linear => Linear
		Procedure, Public ::  RungeKutta=> sch_RungeKutta
		Procedure, Private ::  FRunge=> sch_FRunge
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

		Allocate(this.ku1(f_x: l_x, f_y : l_y))
		Allocate(this.ku2(f_x: l_x, f_y : l_y))
		Allocate(this.ku3(f_x: l_x, f_y : l_y))
		Allocate(this.ku4(f_x: l_x, f_y : l_y))

		Allocate(this.kv1(f_x: l_x, f_y : l_y))
		Allocate(this.kv2(f_x: l_x, f_y : l_y))
		Allocate(this.kv3(f_x: l_x, f_y : l_y))
		Allocate(this.kv4(f_x: l_x, f_y : l_y))

		Allocate(this.kh1(f_x: l_x, f_y : l_y))
		Allocate(this.kh2(f_x: l_x, f_y : l_y))
		Allocate(this.kh3(f_x: l_x, f_y : l_y))
		Allocate(this.kh4(f_x: l_x, f_y : l_y))


End Subroutine



subroutine Linear(this, var, var_pr, grid)

	Class(schema) :: this
	Class(f_var) :: var, var_pr
	Class(g_var) :: grid

	real(8) g, height, dt, partial(1:2), temp(1:3)
	integer(4) face, x, y, dim, i, j, stat, ns_x, ns_y, nf_x, nf_y, ier

	g = grid.g;  height = var_pr.height;  dim = var_pr.dim
	dt = grid.dt;  ns_x = var.ns_x;  ns_y = var.ns_y
	nf_x = var.nf_x;  nf_y = var.nf_y

	do face = 1, 6
		do y = ns_y, nf_y
			do x = ns_x, nf_x

				partial(1) = grid.partial_c1_x(var_pr.h_height(x-1:x+1, y, face), x, y)
				var.x_vel(x, y, face) = var_pr.x_vel(x, y, face) - dt*g*partial(1)

				temp(:) = var_pr.h_height(x, y-1:y+1, face)
				partial(1) = grid.partial_c1_y(temp, x, y)
				var.y_vel(x, y, face) = var_pr.y_vel(x, y, face) - dt*g*partial(1)

				temp(:) = var_pr.x_vel(x-1:x+1, y, face)
				partial(1) = grid.partial_c1_x(temp, x, y)
				temp(:) = var_pr.y_vel(x, y-1:y+1, face)
				partial(2) = grid.partial_c1_y(temp, x, y)
				var.h_height(x, y, face) = var_pr.h_height(x, y, face) - dt*height*(partial(1) + partial(2))


			end do
		end do
	end do

		call var_pr.equal(var)

end subroutine




Subroutine sch_RungeKutta(this, var, var_pr, grid)

	Class(schema) :: this
	Class(f_var) :: var, var_pr
	Class(g_var) :: grid

	real(8) g, height, dt, partial(1:2), temp(1:3)
	integer(4) face, x, y, dim, i, j, stat, ns_x, ns_y, nf_x, nf_y, ier

	g = grid.g;  height = var_pr.height;  dim = var_pr.dim
	dt = grid.dt;  ns_x = var.ns_x;  ns_y = var.ns_y
	nf_x = var.nf_x;  nf_y = var.nf_y


	do face = 1, 6

		call this.FRunge(this.ku1, this.kv1, this.kh1, var_pr, grid, face)
		call this.FRunge(this.ku2, this.kv2, this.kh2, var_pr, grid, face)
		call this.FRunge(this.ku3, this.kv3, this.kh3, var_pr, grid, face)
		call this.FRunge(this.ku4, this.kv4, this.kh4, var_pr, grid, face)

		do y = ns_y, nf_y
			do x= ns_x, nf_x
var.x_vel(x, y, face) = var_pr.x_vel(x, y, face) + (this.ku1(x, y) + 2.0*this.ku2(x, y) + 2.0*this.ku3(x, y) + this.ku4(x, y))/6.0
var.y_vel(x, y, face) = var_pr.y_vel(x, y, face) + (this.kv1(x, y) + 2.0*this.kv2(x, y) + 2.0*this.kv3(x, y) + this.kv4(x, y))/6.0
var.h_height(x, y, face) = var_pr.h_height(x, y, face) + (this.kh1(x, y) + 2.0*this.kh2(x, y) + 2.0*this.kh3(x, y) + this.kh4(x, y))/6.0
			end do
		end do
	end do

	call var_pr.equal(var)

End Subroutine


Subroutine sch_FRunge(this, ku, kv, kh, var, grid, face)
	Class(schema) :: this
	Class(f_var) :: var
	Class(g_var) :: grid

	integer(4), intent(in) :: face
	integer(4) :: x,y

	Real(8), Intent(inout) :: ku(var.first_x : var.last_x, var.first_y : var.last_y)
	Real(8), Intent(inout) :: kv(var.first_x : var.last_x, var.first_y : var.last_y)
	Real(8), Intent(inout) :: kh(var.first_x : var.last_x, var.first_y : var.last_y)

	Real(8) :: g, height, dt, partial(1:2), temp(1:3)

	g = grid.g; height = var.height

	do y = var.ns_y, var.nf_y
		do x = var.ns_x, var.nf_x

			partial(1) = grid.partial_c1_x(var.h_height(x-1:x+1, y, face), x, y)
			ku(x, y) = - dt*g*partial(1)

			temp(:) = var.h_height(x, y-1:y+1, face)
			partial(1) = grid.partial_c1_y(temp, x, y)
			kv(x, y) =  - dt*g*partial(1)

			temp(:) = var.x_vel(x-1:x+1, y, face)
			partial(1) = grid.partial_c1_x(temp, x, y)
			temp(:) = var.y_vel(x, y-1:y+1, face)
			partial(2) = grid.partial_c1_y(temp, x, y)
			kh(x, y) = - dt*height*(partial(1) + partial(2))


		end do
	end do

end Subroutine




Subroutine deinit(this)
	Class(schema) :: this

	if (Allocated(this.ku1)) Deallocate(this.ku1)
	if (Allocated(this.ku2)) Deallocate(this.ku2)
	if (Allocated(this.ku3)) Deallocate(this.ku3)
	if (Allocated(this.ku4)) Deallocate(this.ku4)

	if (Allocated(this.kv1)) Deallocate(this.kv1)
	if (Allocated(this.kv2)) Deallocate(this.kv2)
	if (Allocated(this.kv3)) Deallocate(this.kv3)
	if (Allocated(this.kv4)) Deallocate(this.kv4)

	if (Allocated(this.kh1)) Deallocate(this.kh1)
	if (Allocated(this.kh2)) Deallocate(this.kh2)
	if (Allocated(this.kh3)) Deallocate(this.kh3)
	if (Allocated(this.kh4)) Deallocate(this.kh4)

End Subroutine




end module
