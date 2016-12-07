module methods

	use derivatives, Only: der
	use metrics, Only: metric
	use func_var, Only: f_var
	use interpolation, Only: interp
	use messenger, Only: message
	use mpi
	use omp_lib

	implicit none


	Private
	Public :: method

	Type method

		Integer(4) first_x, first_y, last_x, last_y, step, dim, ns_x, ns_y, nf_x, nf_y
		Real(8) :: g, dt, dh, height

		Real(8), Allocatable :: ku_cov(:, :, :, :)
		Real(8), Allocatable :: kv_cov(:, :, :, :)
		Real(8), Allocatable :: ku_con(:, :, :, :)
		Real(8), Allocatable :: kv_con(:, :, :, :)
		Real(8), Allocatable :: kh(:, :, :, :)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: Euler => Euler
		Procedure, Public :: Predictor_corrector => Predictor_corrector
		Procedure, Public ::  RungeKutta=> RungeKutta
		Procedure, Private ::  FRunge=> FRunge
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS


Subroutine init(this, f)

	Class(f_var) :: f
	Class(method) :: this

	Integer(4) f_x, f_y, l_x, l_y

	f_x = f.first_x;  f_y = f.first_y
	l_x = f.last_x;  l_y = f.last_y

	this.first_x = f_x;  this.first_y = f_y
	this.last_x = l_x;  this.last_y = l_y

	this.step = f.step;  this.dim = f.dim
	this.ns_x = f.ns_x;  this.ns_y = f.ns_y
	this.nf_x = f.nf_x;  this.nf_y = f.nf_y

	this.dt = f.dt;  this.g = f.g;  this.dh = f.delta_on_cube
	this.height = f.height

	Allocate(this.ku_cov(f_x: l_x, f_y : l_y, 6, 0:4))
	Allocate(this.kv_cov(f_x: l_x, f_y : l_y, 6, 0:4))
	Allocate(this.ku_con(f_x: l_x, f_y : l_y, 6, 0:4))
	Allocate(this.kv_con(f_x: l_x, f_y : l_y, 6, 0:4))
	Allocate(this.kh(f_x: l_x, f_y : l_y, 6, 0:4))


End Subroutine



Subroutine Euler(this, var, var_pr, metr, inter, msg)

	Class(method) :: this
	Class(f_var) :: var, var_pr
	Class(metric) :: metr
	Class(interp) :: inter
	Class(message) :: msg
	Type(der) :: d

	Real(8) g, height, dt, partial, temp1(-this.step:this.step), temp2(-this.step:this.step), div, dh
	Integer(4) face, x, y, dim, step

	g = this.g;  height = this.height;  dim = this.dim
	dt = this.dt;  step = this.step;  dh = this.dh

	!$OMP PARALLEL PRIVATE(face, y, x, partial, temp1, temp2, div)
	!$OMP DO

	do face = 1, 6
		do y = var.ns_y, var.nf_y
			do x = var.ns_x, var.nf_x

				temp1(:) = var_pr.h_height(x-step:x+step, y, face)
				partial = d.partial_c(temp1, dh, step)
				var.u_cov(x, y, face) = var_pr.u_cov(x, y, face) - dt*g*partial

				temp1(:) = var_pr.h_height(x, y-step:y+step, face)
				partial = d.partial_c(temp1, dh, step)
				var.v_cov(x, y, face) = var_pr.v_cov(x, y, face) - dt*g*partial

				temp1(:) = var_pr.u_con(x-step:x+step, y, face)
				temp2(:) = var_pr.v_con(x, y-step:y+step, face)
				div = d.div(metr, temp1, temp2, dh, x, y, step)
				var.h_height(x, y, face) = var_pr.h_height(x, y, face) - dt*height*div
			end do
		end do
	end do

	!$OMP END DO
	!$OMP END PARALLEL

	call var_pr.equal(var, metr)
	call msg.msg(var_pr.h_height, var_pr.lon_vel, var_pr.lat_vel)
	call var_pr.interpolate(inter, metr)

end Subroutine



Subroutine Predictor_corrector(this, var, var_pr, metr, inter, msg)

	Class(method) :: this
	Class(f_var) :: var, var_pr
	Class(metric) :: metr
	Class(interp) :: inter
	Class(message) :: msg
	Type(der) :: d

	Real(8) g, height, dt, partial, temp1(-this.step:this.step), temp2(-this.step:this.step), div, dh
	Integer(4) face, x, y, dim, step, i

	g = this.g;  height = this.height;  dim = this.dim;
	dt = this.dt;  step = this.step;  dh = this.dh

	this.ku_con(:,:,:,0) = var_pr.u_con(:,:,:)
	this.kv_con(:,:,:,0) = var_pr.v_con(:,:,:)
	this.ku_cov(:,:,:,0) = var_pr.u_cov(:,:,:)
	this.kv_cov(:,:,:,0) = var_pr.v_cov(:,:,:)
	this.kh(:,:,:,0) = var_pr.h_height(:,:,:)

	call this.Euler(var, var_pr, metr, inter, msg) ! Predictor

	do i = 1, 2 ! Corrector

		!$OMP PARALLEL PRIVATE(face, y, x, partial, temp1, temp2, div)
		!$OMP DO

		do face = 1, 6
			do y = var.ns_y, var.nf_y
				do x = var.ns_x, var.nf_x

					temp1(:) = (var_pr.h_height(x-step:x+step, y, face) + this.kh(x-step:x+step, y, face, 0))/2d0
					partial = d.partial_c(temp1, dh, step)
					var.u_cov(x, y, face) = this.ku_cov(x, y, face, 0) - dt*g*partial

					temp1(:) = (var_pr.h_height(x, y-step:y+step, face) + this.kh(x, y-step:y+step, face, 0))/2d0
					partial = d.partial_c(temp1, dh, step)
					var.v_cov(x, y, face) = this.kv_cov(x, y, face, 0) - dt*g*partial

					temp1(:) = (var_pr.u_con(x-step:x+step, y, face) + this.ku_con(x-step:x+step, y, face, 0))/2d0
					temp2(:) = (var_pr.v_con(x, y-step:y+step, face) + this.kv_con(x, y-step:y+step, face, 0))/2d0
					div = d.div(metr, temp1, temp2, dh, x, y, step)
					var.h_height(x, y, face) = this.kh(x, y, face, 0) - dt*height*div

				end do
			end do
		end do

		!$OMP END DO
		!$OMP END PARALLEL

		call var_pr.equal(var, metr)
		call msg.msg(var_pr.h_height, var_pr.lon_vel, var_pr.lat_vel)
		call var_pr.interpolate(inter, metr)

	end do

end Subroutine




Subroutine RungeKutta(this, var, var_pr, metr, inter, msg)

	Class(method) :: this
	Class(f_var) :: var, var_pr
	Class(metric) :: metr
	Class(interp) :: inter
	Class(message) :: msg

	Integer(4) face, x, y, i, j, stat, ns_x, ns_y, nf_x, nf_y, ier, iteration

	ns_x = this.ns_x;  ns_y = this.ns_y
	nf_x = this.nf_x;  nf_y = this.nf_y

	call var_pr.cov_to_con(metr)


	this.ku_cov(:, :, :, 0) = var_pr.u_cov(:, :, :)
	this.kv_cov(:, :, :, 0) = var_pr.v_cov(:, :, :)
	this.ku_con(:, :, :, 0) = var_pr.u_con(:, :, :)
	this.kv_con(:, :, :, 0) = var_pr.v_con(:, :, :)
	this.kh(:, :, :, 0) = var_pr.h_height(:, :, :)

		do iteration = 1, 4
			call this.FRunge(metr, iteration)
			var_pr.u_cov(:, :, :) = this.ku_cov(:, :, :, iteration)
			var_pr.v_cov(:, :, :) = this.kv_cov(:, :, :, iteration)
			var_pr.h_height(:, :, :) = this.kh(:, :, :, iteration)
			call var_pr.equal(var_pr, metr)
			call msg.msg(var_pr.h_height, var_pr.lon_vel, var_pr.lat_vel)
			call var_pr.interpolate(inter, metr)
			this.ku_cov(:, :, :, iteration) = var_pr.u_cov(:, :, :)
			this.kv_cov(:, :, :, iteration) = var_pr.v_cov(:, :, :)
			this.ku_con(:, :, :, iteration) = var_pr.u_con(:, :, :)
			this.kv_con(:, :, :, iteration) = var_pr.v_con(:, :, :)
			this.kh(:, :, :, iteration) = var_pr.h_height(:, :, :)
		end do

	!$OMP PARALLEL PRIVATE(face, y, x)
	!$OMP DO

	do face = 1, 6
		do y = ns_y, nf_y
			do x= ns_x, nf_x
var.u_cov(x, y, face) = this.ku_cov(x, y, face, 0) + (this.ku_cov(x, y, face, 1) + 2d0*this.ku_cov(x, y, face, 2) + 2d0*this.ku_cov(x, y, face, 3) + this.ku_cov(x, y, face, 4))/6d0
var.v_cov(x, y, face) = this.kv_cov(x, y, face, 0) + (this.kv_cov(x, y, face, 1) + 2d0*this.kv_cov(x, y, face, 2) + 2d0*this.kv_cov(x, y, face, 3) + this.kv_cov(x, y, face, 4))/6d0
var.h_height(x, y, face) = this.kh(x, y, face, 0) + (this.kh(x, y, face, 1) + 2d0*this.kh(x, y, face, 2) + 2d0*this.kh(x, y, face, 3) + this.kh(x, y, face, 4))/6d0
			end do
		end do
	end do

	!$OMP END DO
	!$OMP END PARALLEL


	call var_pr.equal(var, metr)
	call msg.msg(var_pr.h_height, var_pr.lon_vel, var_pr.lat_vel)
	call var_pr.interpolate(inter, metr)

End Subroutine



Subroutine FRunge(this, metr, i)
	Class(method) :: this
	Class(metric) :: metr
	Type(der) :: d

	Integer(4), intent(in) :: i
	Real(8) g, height, dt, partial, temp1(-this.step:this.step), temp2(-this.step:this.step), coef(0:3), div, dh
	Integer(4) x,y, face, step, ns_x, ns_y, nf_x, nf_y

	coef(0) = 0d0;  coef(1) = 5d-1;  coef(2) = 5d-1;  coef(3) = 1d0;

	dt = this.dt;  g = this.g; height = this.height
	dh = this.dh;  step = this.step
	ns_x = this.ns_x;  ns_y = this.ns_y
	nf_x = this.nf_x;  nf_y = this.nf_y

	!$OMP PARALLEL PRIVATE(face, y, x, partial, temp1, temp2, div)
	!$OMP DO

	do face = 1, 6
		do y = ns_y, nf_y
			do x = ns_x, nf_x

				temp1(:) = this.kh(x-step:x+step, y, face, 0) + coef(i-1)*this.kh(x-step:x+step, y, face, i-1)
				partial = d.partial_c(temp1, dh, step)
				this.ku_cov(x, y, face, i) = - dt*g*partial

				temp2(:) = this.kh(x, y-step:y+step, face, 0) + coef(i-1)*this.kh(x, y-step:y+step, face, i-1)
				partial = d.partial_c(temp2, dh, step)
				this.kv_cov(x, y, face, i) =  - dt*g*partial

				temp1(:) = this.ku_con(x-step:x+step, y, face, 0) + coef(i-1)*this.ku_con(x-step:x+step, y, face, i-1)
				temp2(:) = this.kv_con(x, y-step:y+step, face, 0) + coef(i-1)*this.kv_con(x, y-step:y+step, face, i-1)
				div = d.div(metr, temp1, temp2, dh, x, y, step)
				this.kh(x, y, face, i) = - dt*height*div

			end do
		end do
	end do

	!$OMP END DO
	!$OMP END PARALLEL

end Subroutine




Subroutine deinit(this)
	Class(method) :: this

	if (Allocated(this.ku_cov)) Deallocate(this.ku_cov)
	if (Allocated(this.kv_cov)) Deallocate(this.kv_cov)
	if (Allocated(this.ku_con)) Deallocate(this.ku_con)
	if (Allocated(this.kv_con)) Deallocate(this.kv_con)
	if (Allocated(this.kh)) Deallocate(this.kh)

End Subroutine



end module