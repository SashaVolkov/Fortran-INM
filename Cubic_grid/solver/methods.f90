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

		Integer(4) first_x, first_y, last_x, last_y, step, dim, ns_x, ns_y, nf_x, nf_y, test
		Real(8) :: g, dt, dh, height

		Real(8), Allocatable :: ku_con(:, :, :, :)
		Real(8), Allocatable :: kv_con(:, :, :, :)
		Real(8), Allocatable :: kh(:, :, :, :)
		Real(8), Allocatable :: msg_time(:,:)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: Euler => Euler
		Procedure, Public :: Predictor_corrector => Predictor_corrector
		Procedure, Public ::  RungeKutta=> RungeKutta
		Procedure, Private ::  FRunge=> FRunge
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS


Subroutine init(this, f, step, Tmax, dt)

	Class(f_var) :: f
	Class(method) :: this
	Integer(4), intent(in) :: step, Tmax
	Real(8), intent(in) :: dt
	Integer(4) f_x, f_y, l_x, l_y

	f_x = f.first_x;  f_y = f.first_y
	l_x = f.last_x;  l_y = f.last_y

	this.first_x = f_x;  this.first_y = f_y
	this.last_x = l_x;  this.last_y = l_y

	this.step = step;  this.dim = f.dim
	this.ns_x = f.ns_x;  this.ns_y = f.ns_y
	this.nf_x = f.nf_x;  this.nf_y = f.nf_y

	this.dt = dt;  this.g = f.g;  this.dh = f.delta_on_cube
	this.height = f.height;  this.test = f.test

	Allocate(this.ku_con(f_x: l_x, f_y : l_y, 6, 0:4))
	Allocate(this.kv_con(f_x: l_x, f_y : l_y, 6, 0:4))
	Allocate(this.kh(f_x: l_x, f_y : l_y, 6, 0:4))
	Allocate(this.msg_time(5,Tmax))


End Subroutine



Subroutine Euler(this, var, var_pr, metr, inter, msg)

	Class(method) :: this
	Class(f_var) :: var, var_pr
	Class(metric) :: metr
	Class(interp) :: inter
	Class(message) :: msg
	Type(der) :: d

	Real(8) :: temp1(-this.step:this.step, -this.step:this.step), temp2(-this.step:this.step, -this.step:this.step)
	Real(8) :: temp1_cov(-this.step:this.step, -this.step:this.step), temp2_cov(-this.step:this.step, -this.step:this.step)
	Real(8) :: g, height, dt, grad_Fx, grad_Fy, S_c(2), S_p(2), h, u_cov, v_cov, div, dh, uu(2)
	Integer(4) x,y, face, step, ns_x, ns_y, nf_x, nf_y

	g = this.g;  height = this.height;
	dt = this.dt;  step = this.step;  dh = this.dh

	!$OMP PARALLEL PRIVATE(face, y, x, grad_Fy, grad_Fx, temp1, temp2, temp1_cov, temp2_cov, div, h, u_cov, v_cov, uu, S_p, S_c)
	!$OMP DO

	do face = 1, 6
		do y = var.ns_y, var.nf_y
			do x = var.ns_x, var.nf_x

			uu = 0d0
			temp1 = var_pr.h_height(x-step:x+step, y-step:y+step, face)
			h = temp1(0,0)
			grad_Fx = d.partial_c(temp1(:,0), dh, step)
			grad_Fy = d.partial_c(temp1(0,:), dh, step)

			temp1 = var_pr.u_con(x-step:x+step, y-step:y+step, face)
			temp2 = var_pr.v_con(x-step:x+step, y-step:y+step, face)

			uu(1) = (temp1(0,0)**2)*metr.Christoffel_x1(1,1,x,y) + 2d0*temp1(0,0)*temp2(0,0)*metr.Christoffel_x1(2,1,x,y)
			uu(2) = (temp2(0,0)**2)*metr.Christoffel_x2(2,2,x,y) + 2d0*temp1(0,0)*temp2(0,0)*metr.Christoffel_x2(2,1,x,y)
			! print *, uu, x, y, face

			uu(1) = dt*(temp1(0,0)*d.partial_c(temp1(:,0), dh, step) + temp2(0,0)*d.partial_c(temp1(0,:), dh, step) + uu(1))
			uu(2) = dt*(temp1(0,0)*d.partial_c(temp2(:,0), dh, step) + temp2(0,0)*d.partial_c(temp2(0,:), dh, step) + uu(2))


			div = d.div(metr, temp1(:,0), temp2(0,:), dh, x, y, step)

			S_c(1) = dt*metr.G_sqr(x, y)*temp2(0,0)*(var.f(x, y, face))
			S_c(2) = dt*metr.G_sqr(x, y)*temp1(0,0)*(var.f(x, y, face))

			S_p(1) = - dt*g*grad_Fx
			S_p(2) = - dt*g*grad_Fy
			call metr.cov_to_con(S_p(1), S_p(2), S_p(1), S_p(2), x, y)
			call metr.cov_to_con(S_c(1), - S_c(2), S_c(1), S_c(2), x, y)

			var.u_con(x, y, face) = var_pr.u_con(x, y, face) - uu(1) + S_c(1) + S_p(1)
			var.v_con(x, y, face) = var_pr.v_con(x, y, face) - uu(2) + S_c(2) + S_p(2)
			var.h_height(x, y, face) = var_pr.h_height(x, y, face) - dt*h*div - temp1(0,0)*dt*grad_Fx - temp2(0,0)*dt*grad_Fy

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

	Real(8) :: temp1(-this.step:this.step, -this.step:this.step), temp2(-this.step:this.step, -this.step:this.step)
	Real(8) :: temp1_cov(-this.step:this.step, -this.step:this.step), temp2_cov(-this.step:this.step, -this.step:this.step)
	Real(8) :: g, height, dt, grad_Fx, grad_Fy, S_c(2), S_p(2), h, u_cov, v_cov, div, dh, uu(2)
	Integer(4) x,y, i, face, step, ns_x, ns_y, nf_x, nf_y
	g = this.g;  height = this.height;
	dt = this.dt;  step = this.step;  dh = this.dh

	this.ku_con(:,:,:,0) = var_pr.u_con(:,:,:)
	this.kv_con(:,:,:,0) = var_pr.v_con(:,:,:)
	this.kh(:,:,:,0) = var_pr.h_height(:,:,:)

	call this.Euler(var, var_pr, metr, inter, msg) ! Predictor

	do i = 1, 2 ! Corrector

		!$OMP PARALLEL PRIVATE(face, y, x, grad_Fy, grad_Fx, temp1, temp2, temp1_cov, temp2_cov, div, h, u_cov, v_cov, uu, S_p, S_c)
		!$OMP DO

		do face = 1, 6
			do y = var.ns_y, var.nf_y
				do x = var.ns_x, var.nf_x

				uu = 0d0
				temp1 = (var_pr.h_height(x-step:x+step, y-step:y+step, face) + this.kh(x, y, face, 0))/2d0
				h = temp1(0,0)
				grad_Fx = d.partial_c(temp1(:,0), dh, step)
				grad_Fy = d.partial_c(temp1(0,:), dh, step)

				temp1 = (var_pr.u_con(x-step:x+step, y-step:y+step, face) + this.ku_con(x, y, face, 0))/2d0
				temp2 = (var_pr.v_con(x-step:x+step, y-step:y+step, face) + this.kv_con(x, y, face, 0))/2d0

				uu(1) = (temp1(0,0)**2)*metr.Christoffel_x1(1,1,x,y) + 2d0*temp1(0,0)*temp2(0,0)*metr.Christoffel_x1(2,1,x,y)
				uu(2) = (temp2(0,0)**2)*metr.Christoffel_x2(2,2,x,y) + 2d0*temp1(0,0)*temp2(0,0)*metr.Christoffel_x2(2,1,x,y)
				! print *, uu, x, y, face

				uu(1) = dt*(temp1(0,0)*d.partial_c(temp1(:,0), dh, step) + temp2(0,0)*d.partial_c(temp1(0,:), dh, step) + uu(1))
				uu(2) = dt*(temp1(0,0)*d.partial_c(temp2(:,0), dh, step) + temp2(0,0)*d.partial_c(temp2(0,:), dh, step) + uu(2))

				div = d.div(metr, temp1(:,0), temp2(0,:), dh, x, y, step)

				S_c(1) = dt*metr.G_sqr(x, y)*temp2(0,0)*(var.f(x, y, face))
				S_c(2) = dt*metr.G_sqr(x, y)*temp1(0,0)*(var.f(x, y, face))

				S_p(1) = - dt*g*grad_Fx
				S_p(2) = - dt*g*grad_Fy
				call metr.cov_to_con(S_p(1), S_p(2), S_p(1), S_p(2), x, y)
				call metr.cov_to_con(S_c(1), - S_c(2), S_c(1), S_c(2), x, y)

				var.u_con(x, y, face) = this.ku_con(x, y, face, 0) - uu(1) + S_c(1) + S_p(1)
				var.v_con(x, y, face) = this.kv_con(x, y, face, 0) - uu(2) + S_c(2) + S_p(2)
				var.h_height(x, y, face) = this.kh(x, y, face, 0) - dt*h*div - temp1(0,0)*dt*grad_Fx - temp2(0,0)*dt*grad_Fy

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




Subroutine RungeKutta(this, var, var_pr, metr, inter, msg, time)

	Class(method) :: this
	Class(f_var) :: var, var_pr
	Class(metric) :: metr
	Class(interp) :: inter
	Class(message) :: msg
	Integer(4), intent(in) :: time

	Integer(4) face, x, y, i, j, stat, ns_x, ns_y, nf_x, nf_y, ier, iteration

	ns_x = this.ns_x;  ns_y = this.ns_y
	nf_x = this.nf_x;  nf_y = this.nf_y


	this.ku_con(:, :, :, 0) = var_pr.u_con(:, :, :)
	this.kv_con(:, :, :, 0) = var_pr.v_con(:, :, :)
	this.kh(:, :, :, 0) = var_pr.h_height(:, :, :)

		do iteration = 1, 4
			call this.FRunge(metr, var_pr, iteration)
			var_pr.u_con(:, :, :) = this.ku_con(:, :, :, iteration)
			var_pr.v_con(:, :, :) = this.kv_con(:, :, :, iteration)
			var_pr.h_height(:, :, :) = this.kh(:, :, :, iteration)
			call var_pr.equal(var_pr, metr)
			call MPI_Barrier(MPI_COMM_WORLD, ier)
			this.msg_time(iteration,time) = MPI_Wtime()
			call msg.msg(var_pr.h_height, var_pr.lon_vel, var_pr.lat_vel)
			call MPI_Barrier(MPI_COMM_WORLD, ier)
			this.msg_time(iteration,time) = MPI_Wtime() - this.msg_time(iteration,time)
			call var_pr.interpolate(inter, metr)
			this.ku_con(:, :, :, iteration) = var_pr.u_con(:, :, :)
			this.kv_con(:, :, :, iteration) = var_pr.v_con(:, :, :)
			this.kh(:, :, :, iteration) = var_pr.h_height(:, :, :)
		end do


	do face = 1, 6
		do y = ns_y, nf_y
			do x= ns_x, nf_x
var.u_con(x, y, face) = this.ku_con(x, y, face, 0) + (this.ku_con(x, y, face, 1) + 2d0*this.ku_con(x, y, face, 2) + 2d0*this.ku_con(x, y, face, 3) + this.ku_con(x, y, face, 4))/6d0
var.v_con(x, y, face) = this.kv_con(x, y, face, 0) + (this.kv_con(x, y, face, 1) + 2d0*this.kv_con(x, y, face, 2) + 2d0*this.kv_con(x, y, face, 3) + this.kv_con(x, y, face, 4))/6d0
var.h_height(x, y, face) = this.kh(x, y, face, 0) + (this.kh(x, y, face, 1) + 2d0*this.kh(x, y, face, 2) + 2d0*this.kh(x, y, face, 3) + this.kh(x, y, face, 4))/6d0
			end do
		end do
	end do


	call var_pr.equal(var, metr)
	call MPI_Barrier(MPI_COMM_WORLD, ier)
	this.msg_time(5,time) = MPI_Wtime()
	call msg.msg(var_pr.h_height, var_pr.lon_vel, var_pr.lat_vel)
	call MPI_Barrier(MPI_COMM_WORLD, ier)
	this.msg_time(5,time) = MPI_Wtime() - this.msg_time(5,time)
	call var_pr.interpolate(inter, metr)

End Subroutine



Subroutine FRunge(this, metr, var, i)
	Class(method) :: this
	Class(metric) :: metr
	Class(f_var) :: var
	Type(der) :: d

	Integer(4), intent(in) :: i
	Real(8) :: temp1(-this.step:this.step, -this.step:this.step), temp2(-this.step:this.step, -this.step:this.step)
	Real(8) :: temp1_cov(-this.step:this.step, -this.step:this.step), temp2_cov(-this.step:this.step, -this.step:this.step)
	Real(8) :: g, height, dt, grad_Fx, grad_Fy, S_c(2), S_p(2), h, u_cov, v_cov, coef(0:3), div, dh, uu(2)
	Integer(4) x,y, face, step, ns_x, ns_y, nf_x, nf_y

	coef(0) = 0d0;  coef(1) = 5d-1;  coef(2) = 5d-1;  coef(3) = 1d0;

	dt = this.dt;  g = this.g; height = this.height
	dh = this.dh;  step = this.step
	ns_x = this.ns_x;  ns_y = this.ns_y
	nf_x = this.nf_x;  nf_y = this.nf_y

	!$OMP PARALLEL PRIVATE(face, y, x, grad_Fy, grad_Fx, temp1, temp2, temp1_cov, temp2_cov, div, h, u_cov, v_cov, uu, S_p, S_c)
	!$OMP DO

	do face = 1,6
		do y = ns_y, nf_y
			do x = ns_x, nf_x

uu = 0d0
temp1 = this.kh(x-step:x+step, y-step:y+step, face, 0) + coef(i-1)*this.kh(x-step:x+step, y-step:y+step, face, i-1)
h = temp1(0,0)
! height = var.h_depth(x, y, face)
grad_Fx = d.partial_c(temp1(:,0), dh, step)
grad_Fy = d.partial_c(temp1(0,:), dh, step)


temp1 = this.ku_con(x-step:x+step, y-step:y+step, face, 0) + coef(i-1)*this.ku_con(x-step:x+step, y-step:y+step, face, i-1)
temp2 = this.kv_con(x-step:x+step, y-step:y+step, face, 0) + coef(i-1)*this.kv_con(x-step:x+step, y-step:y+step, face, i-1)
if (this.test /= 1) then
uu(1) = (temp1(0,0)**2)*metr.Christoffel_x1(1,1,x,y) + 2d0*temp1(0,0)*temp2(0,0)*metr.Christoffel_x1(2,1,x,y)
uu(2) = (temp2(0,0)**2)*metr.Christoffel_x2(2,2,x,y) + 2d0*temp1(0,0)*temp2(0,0)*metr.Christoffel_x2(2,1,x,y)

uu(1) = dt*(temp1(0,0)*d.partial_c(temp1(:,0), dh, step) + temp2(0,0)*d.partial_c(temp1(0,:), dh, step) + uu(1))
uu(2) = dt*(temp1(0,0)*d.partial_c(temp2(:,0), dh, step) + temp2(0,0)*d.partial_c(temp2(0,:), dh, step) + uu(2))

div = d.div(metr, temp1(:,0), temp2(0,:), dh, x, y, step)

S_c(1) = dt*metr.G_sqr(x, y)*temp2(0,0)*(var.f(x, y, face))
S_c(2) = dt*metr.G_sqr(x, y)*temp1(0,0)*(var.f(x, y, face))

S_p(1) = - dt*g*grad_Fx
S_p(2) = - dt*g*grad_Fy
call metr.cov_to_con(S_p(1), S_p(2), S_p(1), S_p(2), x, y)
call metr.cov_to_con(S_c(1), - S_c(2), S_c(1), S_c(2), x, y)

this.ku_con(x, y, face, i) = - uu(1) + S_c(1) + S_p(1)
this.kv_con(x, y, face, i) = - uu(2) + S_c(2) + S_p(2)
end if
this.kh(x, y, face, i) = - dt*h*div - temp1(0,0)*dt*grad_Fx - temp2(0,0)*dt*grad_Fy


			end do
		end do
	end do

	!$OMP END DO
	!$OMP END PARALLEL

end Subroutine




Subroutine deinit(this)
	Class(method) :: this

	if (Allocated(this.ku_con)) Deallocate(this.ku_con)
	if (Allocated(this.kv_con)) Deallocate(this.kv_con)
	if (Allocated(this.kh)) Deallocate(this.kh)
	if (Allocated(this.msg_time)) Deallocate(this.msg_time)

End Subroutine



end module