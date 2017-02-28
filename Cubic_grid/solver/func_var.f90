module func_var

	use sphere_geometry, Only: geometry
	use interpolation, Only: interp
	use metrics, Only: metric

implicit none

	Private
	Public :: f_var

	Type f_var

		Real(8), Allocatable :: h_height(:, :, :)
		Real(8), Allocatable :: h_starter(:, :, :)
		Real(8), Allocatable :: u_cov(:, :, :)
		Real(8), Allocatable :: v_cov(:, :, :)
		Real(8), Allocatable :: u_con(:, :, :)
		Real(8), Allocatable :: v_con(:, :, :)
		Real(8), Allocatable :: lon_vel(:, :, :)
		Real(8), Allocatable :: lat_vel(:, :, :)
		Real(8), Allocatable :: f_cor(:, :, :)
		Real(8) height,  g, dt, delta_on_cube, u0
		Integer(4) step, dim, interp_factor(1:4), Neighbours_face(6, 4), grid_type, rescale
		Integer(4) ns_x, ns_y, nf_x, nf_y, first_x, first_y, last_x, last_y, snd_xy(6, 4, 2), rcv_xy(6, 4, 2)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Public :: equal => equal
		Procedure, Public :: start_conditions => start_conditions
		Procedure, Public :: interpolate => interpolate
		Procedure, Private :: Velocity_from_spherical => Velocity_from_spherical
		Procedure, Public :: Velocity_to_spherical => Velocity_to_spherical
		Procedure, Private :: Velocity_from_spherical_border => Velocity_from_spherical_border
		Procedure, Private :: Velocity_to_spherical_border => Velocity_to_spherical_border
		Procedure, Public :: cov_to_con => cov_to_con
		Procedure, Public :: con_to_cov => con_to_cov
	End Type


CONTAINS



	Subroutine init(this, metr, height, dt)

		Class(f_var) :: this
		Class(metric) :: metr
		Real(8), intent(in) :: height, dt
		Integer(4) :: i, x, y, face

		this.g = 980616d-5;  this.dt = dt

		this.ns_x = metr.ns_xy(1);  this.ns_y = metr.ns_xy(2)
		this.nf_x = metr.nf_xy(1);  this.nf_y = metr.nf_xy(2)

		this.first_x = metr.first_x;  this.first_y = metr.first_y
		this.last_x = metr.last_x;  this.last_y = metr.last_y

		this.step = metr.step;  this.height = 294d2/this.g;  this.dim = metr.dim
		this.Neighbours_face = metr.Neighbours_face;  this.grid_type = metr.grid_type

		this.snd_xy = metr.snd_xy;  this.rcv_xy = metr.rcv_xy
		this.interp_factor(:) = 0;  this.rescale = metr.rescale
		this.delta_on_cube = metr.delta_on_cube


		do i = 1, 4
			if(this.grid_type == 1) then
				if(this.Neighbours_face(2, i) /= 2) this.interp_factor(i) = 1
			end if
		end do

		call this.alloc()

	end Subroutine



	Subroutine alloc(this)

		Class(f_var) :: this

		Allocate(this.h_height(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.h_starter(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.u_cov(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.v_cov(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.u_con(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.v_con(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.lon_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.lat_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.f_cor(this.first_x:this.last_x, this.first_y:this.last_y, 6))

	end Subroutine



	Subroutine deinit(this)
		Class(f_var) :: this

		if (Allocated(this.h_height)) Deallocate(this.h_height)
		if (Allocated(this.h_starter)) Deallocate(this.h_starter)
		if (Allocated(this.u_cov)) Deallocate(this.u_cov)
		if (Allocated(this.v_cov)) Deallocate(this.v_cov)
		if (Allocated(this.u_con)) Deallocate(this.u_con)
		if (Allocated(this.v_con)) Deallocate(this.v_con)
		if (Allocated(this.lon_vel)) Deallocate(this.lon_vel)
		if (Allocated(this.lat_vel)) Deallocate(this.lat_vel)
		if (Allocated(this.f_cor)) Deallocate(this.f_cor)

	end Subroutine



	Subroutine equal(var_pr, var, metr)

		Class(f_var) :: var_pr, var
		Class(metric) :: metr


		var_pr.h_height = var.h_height
		var_pr.u_cov = var.u_cov
		var_pr.v_cov = var.v_cov

		call var_pr.Velocity_to_spherical_border(metr)


	end Subroutine



	Subroutine start_conditions(this, metr, geom, omega_cor)

		Class(f_var) :: this
		Class(metric) :: metr
		Class(geometry) :: geom
		Real(8), intent(in) :: omega_cor
		Integer(4) dim, x, y, face
		Real(8) gh0, r, R_BIG, zero(2), pi, u0, alpha, lon, lat

		pi = 314159265358979323846d-20

		dim = this.dim; R_BIG = geom.radius/3d0
		zero(:) = (/0d0, 0d-1*pi/)
		this.h_height = 0d0
		this.u_cov = 0d0;  this.v_cov = 0d0;  this.lon_vel = 0d0
		this.u_con = 0d0;  this.v_con = 0d0;  this.lat_vel = 0d0

		u0 = 2d0*pi*geom.radius/(12d0*24d0*60d0*60d0);  alpha = 0d0 !-pi/4d0
		gh0 = 294d2

		do face = 1, 6
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x

lat = metr.latlon_c(1,x,y,face)
lon = metr.latlon_c(2,x,y,face)

this.lon_vel(x, y, face) = u0
this.lat_vel(x, y, face) = 0d0
this.h_height(x, y, face) = geom.radius*2d0*omega_cor*u0 *(dcos(lat))/this.g
this.f_cor(x, y, face)= (2d0*omega_cor)*dsin(lat)


! r = geom.dist(zero(:),metr.latlon_c(1:2,x,y,face))
! this.h_height(x, y, face) = 10.0*exp(-((r/R_BIG)**2)*4.0)

! this.f_cor(x, y, face)= -2*omega_cor*dsin(lat) ! function of latitude
! if ( face == 2 .and. x == dim ) print *, this.f_cor(x, y, face)
				end do
			end do
		end do

		call this.Velocity_from_spherical(metr)
		call this.con_to_cov(metr)
		this.h_starter(:,:,:) = this.h_height(:,:,:)

	end Subroutine




	Subroutine interpolate(this, i, metr)
		Class(f_var) :: this
		Class(interp) :: i
		Class(metric) :: metr

		call i.Lagrange(this.h_height, this.interp_factor)
		call i.Lagrange(this.lat_vel, this.interp_factor)
		call i.Lagrange(this.lon_vel, this.interp_factor)
		call this.Velocity_from_spherical_border(metr)
! 		this.lat_vel = 0d0
! 		call this.Velocity_from_spherical(metr)
		call this.cov_to_con(metr)

	end Subroutine



	Subroutine cov_to_con(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Real(8) :: G_inv(2,2)
		Integer(4) :: x, y, face

		do face = 1, 6
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x

G_inv(:,:) = metr.G_inverse(:,:, x, y)
this.u_con(x, y, face) = G_inv(1, 1) * this.u_cov(x, y, face) + G_inv(1, 2) * this.v_cov(x, y, face)
this.v_con(x, y, face) = G_inv(2, 2) * this.v_cov(x, y, face) + G_inv(2, 1) * this.u_cov(x, y, face)

				end do
			end do
		end do

	end Subroutine



	Subroutine con_to_cov(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Real(8) :: G(2,2)
		Integer(4) :: x, y, face

		do face = 1, 6
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x

G(:,:) = metr.G_tensor(:,:, x, y)
this.u_cov(x, y, face) = G(1, 1) * this.u_con(x, y, face) + G(1, 2) * this.v_con(x, y, face)
this.v_cov(x, y, face) = G(2, 2) * this.v_con(x, y, face) + G(2, 1) * this.u_con(x, y, face)

				end do
			end do
		end do

	end Subroutine



	Subroutine Velocity_to_spherical(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Real(8) :: vel_x_contr, vel_y_contr, A(2,2)
		Integer(4) :: x, y, face

		do face = 1, 6
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x

A(:,:) = metr.Tr_to_sph(:,:, x, y, face)
this.lon_vel(x, y, face) = A(1, 1) * this.u_con(x, y, face) + A(1, 2) * this.v_con(x, y, face)
this.lat_vel(x, y, face) = A(2, 2) * this.v_con(x, y, face) + A(2, 1) * this.u_con(x, y, face)

				end do
			end do
		end do

	end Subroutine



	Subroutine Velocity_from_spherical(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Real(8) :: A_inv(2,2)
		Integer(4) :: x, y, face

		do face = 1, 6
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x

A_inv(:,:) = metr.Tr_to_cube(:,:, x, y, face)
this.u_con(x, y, face) = A_inv(1, 1) * this.lon_vel(x, y, face) + A_inv(1, 2) * this.lat_vel(x, y, face)
this.v_con(x, y, face) = A_inv(2, 2) * this.lat_vel(x, y, face) + A_inv(2, 1) * this.lon_vel(x, y, face)

				end do
			end do
		end do

	end Subroutine



	Subroutine Velocity_to_spherical_border(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Real(8) :: G_inv(2,2), A(2,2)
		Integer(4) :: x, y, face, i, x_fin(4), y_fin(4)

		x_fin(:) = this.nf_x;  x_fin(4) = this.snd_xy(2, 4, 1) + this.step - 1
		y_fin(:) = this.nf_y;  y_fin(3) = this.snd_xy(2, 3, 2) + this.step - 1

		do face = 1, 6
			do i = 1, 4
				do y = this.snd_xy(2, i, 2), y_fin(i)
					do x = this.snd_xy(2, i, 1), x_fin(i)

					A(:,:) = metr.Tr_to_sph(:,:, x, y, face)
					G_inv(:,:) = metr.G_inverse(:,:, x, y)

this.u_con(x, y, face) = G_inv(1, 1) * this.u_cov(x, y, face) + G_inv(1, 2) * this.v_cov(x, y, face)
this.v_con(x, y, face) = G_inv(2, 2) * this.v_cov(x, y, face) + G_inv(2, 1) * this.u_cov(x, y, face)

this.lon_vel(x, y, face) = A(1, 1) * this.u_con(x, y, face) + A(1, 2) * this.v_con(x, y, face)
this.lat_vel(x, y, face) = A(2, 2) * this.v_con(x, y, face) + A(2, 1) * this.u_con(x, y, face)

					end do
				end do
			end do
		end do

	end Subroutine



	Subroutine Velocity_from_spherical_border(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Real(8) :: G(2,2), A_inv(2,2)
		Integer(4) :: x, y, face, i, x_fin(4), y_fin(4)

		x_fin(:) = this.nf_x;  x_fin(4) = this.rcv_xy(2, 4, 1) + this.step - 1;  x_fin(2) = this.last_x
		y_fin(:) = this.nf_y;  y_fin(3) = this.rcv_xy(2, 3, 2) + this.step - 1;  y_fin(1) = this.last_y

		do face = 1, 6
			do i = 1, 4
				do y = this.rcv_xy(2, i, 2), y_fin(i)
					do x = this.rcv_xy(2, i, 1), x_fin(i)

					A_inv(:,:) = metr.Tr_to_cube(:,:, x, y, face)
					G(:,:) = metr.G_tensor(:,:, x, y)

this.u_con(x, y, face) = A_inv(1, 1) * this.lon_vel(x, y, face) + A_inv(1, 2) * this.lat_vel(x, y, face)
this.v_con(x, y, face) = A_inv(2, 2) * this.lat_vel(x, y, face) + A_inv(2, 1) * this.lon_vel(x, y, face)

this.u_cov(x, y, face) = G(1, 1) * this.u_con(x, y, face) + G(1, 2) * this.v_con(x, y, face)
this.v_cov(x, y, face) = G(2, 2) * this.v_con(x, y, face) + G(2, 1) * this.u_con(x, y, face)

					end do
				end do
			end do
		end do

	end Subroutine




end module