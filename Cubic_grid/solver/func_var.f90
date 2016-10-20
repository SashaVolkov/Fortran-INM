module func_var

	use parallel_cubic, Only: parallel
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
		real(8) height
		integer(4) step, dim, Xsize, Ysize, interp_factor(1:4), Neighbours_face(6, 4)
		integer(4) ns_x, ns_y, nf_x, nf_y, first_x, first_y, last_x, last_y, snd_xy(6, 4, 2), rcv_xy(6, 4, 2)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Public :: equal => equal
		Procedure, Public :: start_conditions => start_conditions
		Procedure, Public :: interpolate => interpolate
		Procedure, Public :: Velocity_from_spherical => Velocity_from_spherical
		Procedure, Public :: Velocity_to_spherical => Velocity_to_spherical
		Procedure, Public :: Velocity_from_spherical_border => Velocity_from_spherical_border
		Procedure, Public :: Velocity_to_spherical_border => Velocity_to_spherical_border
		Procedure, Public :: cov_to_con => cov_to_con
		Procedure, Public :: con_to_cov => con_to_cov
	End Type


CONTAINS



	subroutine init(this, paral, metr, height)

		Class(f_var) :: this
		Class(metric) :: metr
		Class(parallel) :: paral
		real(8), intent(in) :: height
		integer(4) :: i

		this.ns_x = paral.ns_xy(1);  this.ns_y = paral.ns_xy(2)
		this.nf_x = paral.nf_xy(1);  this.nf_y = paral.nf_xy(2)

		this.first_x = paral.first_x;  this.first_y = paral.first_y
		this.last_x = paral.last_x;  this.last_y = paral.last_y

		this.Xsize = paral.Xsize;  this.Ysize = paral.Ysize
		this.step = paral.step;  this.height = height;  this.dim = paral.dim

		this.Neighbours_face = paral.Neighbours_face

		this.snd_xy(:,:,:) = paral.snd_xy(:,:,:)
		this.rcv_xy(:,:,:) = paral.rcv_xy(:,:,:)

		this.interp_factor(:) = 0

		do i = 1, 4
			if(metr.grid_type == 1) then
				if(paral.Neighbours_face(2, i) /= 2) this.interp_factor(i) = 1
			else
				this.interp_factor(i) = 0
			end if
		end do

		call this.alloc()

	end subroutine



	subroutine alloc(this)

		Class(f_var) :: this

		Allocate(this.h_height(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.h_starter(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.u_cov(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.v_cov(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.u_con(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.v_con(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.lon_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.lat_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6))

	end subroutine



	subroutine deinit(this)
		Class(f_var) :: this

		if (Allocated(this.h_height)) Deallocate(this.h_height)
		if (Allocated(this.h_starter)) Deallocate(this.h_starter)
		if (Allocated(this.u_cov)) Deallocate(this.u_cov)
		if (Allocated(this.v_cov)) Deallocate(this.v_cov)
		if (Allocated(this.u_con)) Deallocate(this.u_con)
		if (Allocated(this.v_con)) Deallocate(this.v_con)
		if (Allocated(this.lon_vel)) Deallocate(this.lon_vel)
		if (Allocated(this.lat_vel)) Deallocate(this.lat_vel)

	end subroutine



	subroutine equal(var_pr, var, metr, vec_only)

		Class(f_var) :: var_pr, var
		Class(metric) :: metr

		integer(4), intent(in) :: vec_only

				if(vec_only == 0) var_pr.h_height(:, :, :)=var.h_height(:, :, :)
				if(vec_only == 1) then
					var_pr.u_cov(:, :, :)=var.u_cov(:, :, :)
					var_pr.v_cov(:, :, :)=var.v_cov(:, :, :)

					if(metr.grid_type == 1) then
						call var_pr.Velocity_to_spherical_border(metr)
					else
						var_pr.lon_vel(:, :, :)=var_pr.u_cov(:, :, :)
						var_pr.lat_vel(:, :, :)=var_pr.v_cov(:, :, :)
					end if
				end if


	end subroutine



	subroutine start_conditions(this, metr, geom)

		Class(f_var) :: this
		Class(metric) :: metr
		Class(geometry) :: geom
		integer(4) dim
		real(8) h0, r, R_BIG

		integer(4) x, y, face

		h0 = this.height;  dim = this.dim; R_BIG = geom.radius/3d0

		do face = 1, 6

			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					this.h_height(x, y, face) = 0d0
					this.u_cov(x, y, face) = 0d0
					this.v_cov(x, y, face) = 0d0
				end do
			end do

			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					r = geom.dist((/0d0,0d0/),metr.latlon_c(:,x,y,face))
					if ( r < R_BIG ) then
						this.h_height(x, y, face) = (h0*0.5)*(1d0 + dcos(geom.pi * r/R_BIG))
					else
						this.h_height(x, y, face) = 0d0
					end if
				end do
			end do


			! do y = this.first_y, this.last_y
			! 	do x = this.first_x, this.last_x
			! 		r = geom.dist((/0d0,0d0/),metr.latlon_c(:,x,y,face))
			! 		this.h_height(x, y, face) = h0*exp(-((10d0*r/geom.radius)**2))
			! 	end do
			! end do

! 			do y = this.first_y, this.last_y
! 				do x = this.first_x, this.last_x
! 					this.lon_vel(x,y,face) = 100d0
! 				end do
! 			end do

		end do

! 		call this.Velocity_from_spherical(metr)
! 		call this.con_to_cov(metr)

		this.h_starter(:,:,:) = this.h_height(:,:,:)

	end subroutine




	subroutine interpolate(this, i, metr, vec_only)
		Class(f_var) :: this
		Class(interp) :: i
		Class(metric) :: metr
		integer(4), intent(in) :: vec_only
		if(metr.grid_type == 1) then
			if(vec_only == 0) then
				call i.Lagrange(this.h_height, this.interp_factor)
			else
				call i.Lagrange(this.lat_vel, this.interp_factor)
				call i.Lagrange(this.lon_vel, this.interp_factor)
				call this.Velocity_from_spherical_border(metr)
			end if
		else if (metr.grid_type == 0) then
			if(vec_only == 1) then
				this.u_cov(:, :, :)=this.lon_vel(:, :, :)
				this.v_cov(:, :, :)=this.lat_vel(:, :, :)
			end if
		end if

		if(vec_only == 1) call this.cov_to_con(metr)

	end subroutine



	subroutine cov_to_con(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Integer(4) :: x, y, face

		do face = 1, 6
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x

this.u_con(x, y, face) = metr.G_inverse(1, 1, x, y) * this.u_cov(x, y, face) + metr.G_inverse(1, 2, x, y) * this.v_cov(x, y, face)
this.v_con(x, y, face) = metr.G_inverse(2, 2, x, y) * this.v_cov(x, y, face) + metr.G_inverse(2, 1, x, y) * this.u_cov(x, y, face)

				end do
			end do
		end do

	end subroutine



	subroutine con_to_cov(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Integer(4) :: x, y, face

		do face = 1, 6
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x

this.u_cov(x, y, face) = metr.G_tensor(1, 1, x, y) * this.u_con(x, y, face) + metr.G_tensor(1, 2, x, y) * this.v_con(x, y, face)
this.v_cov(x, y, face) = metr.G_tensor(2, 2, x, y) * this.v_con(x, y, face) + metr.G_tensor(2, 1, x, y) * this.u_con(x, y, face)

				end do
			end do
		end do

	end subroutine



	subroutine Velocity_to_spherical(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Real(8) :: vel_x_contr, vel_y_contr
		Integer(4) :: x, y, face

		do face = 1, 6
			do y = this.ns_y, this.nf_y
				do x = this.ns_x, this.nf_x

this.lon_vel(x, y, face) = metr.J_to_sph(1, 1, x, y, face) * this.u_con(x, y, face) + metr.J_to_sph(1, 2, x, y, face) * this.v_con(x, y, face)
this.lat_vel(x, y, face) = metr.J_to_sph(2, 2, x, y, face) * this.v_con(x, y, face) + metr.J_to_sph(2, 1, x, y, face) * this.u_con(x, y, face)

				end do
			end do
		end do

	end subroutine



	subroutine Velocity_from_spherical(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Integer(4) :: x, y, face

		do face = 1, 6
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x

this.u_con(x, y, face) = metr.J_to_cube(1, 1, x, y, face) * this.lon_vel(x, y, face) + metr.J_to_cube(1, 2, x, y, face) * this.lat_vel(x, y, face)
this.v_con(x, y, face) = metr.J_to_cube(2, 2, x, y, face) * this.lat_vel(x, y, face) + metr.J_to_cube(2, 1, x, y, face) * this.lon_vel(x, y, face)

				end do
			end do
		end do

	end subroutine



	subroutine Velocity_to_spherical_border(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Real(8) :: vel_x_contr, vel_y_contr
		Integer(4) :: x, y, face, i, x_fin(4), y_fin(4)

		x_fin(:) = this.nf_x;  x_fin(4) = this.snd_xy(2, 4, 1) + this.step - 1
		y_fin(:) = this.nf_y;  y_fin(3) = this.snd_xy(2, 3, 2) + this.step - 1

		do face = 1, 6
			do i = 1, 4
				! if(this.interp_factor(i) == 1) then
					do y = this.snd_xy(2, i, 2), y_fin(i)
						do x = this.snd_xy(2, i, 1), x_fin(i)

this.u_con(x, y, face) = metr.G_inverse(1, 1, x, y) * this.u_cov(x, y, face) + metr.G_inverse(1, 2, x, y) * this.v_cov(x, y, face)
this.v_con(x, y, face) = metr.G_inverse(2, 2, x, y) * this.v_cov(x, y, face) + metr.G_inverse(2, 1, x, y) * this.u_cov(x, y, face)

this.lon_vel(x, y, face) = metr.J_to_sph(1, 1, x, y, face) * this.u_con(x, y, face) + metr.J_to_sph(1, 2, x, y, face) * this.v_con(x, y, face)
this.lat_vel(x, y, face) = metr.J_to_sph(2, 2, x, y, face) * this.v_con(x, y, face) + metr.J_to_sph(2, 1, x, y, face) * this.u_con(x, y, face)

						end do
					end do
				! end if
			end do
		end do

	end subroutine



	subroutine Velocity_from_spherical_border(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Integer(4) :: x, y, face, i, x_fin(4), y_fin(4)

		x_fin(:) = this.nf_x;  x_fin(4) = this.rcv_xy(2, 4, 1) + this.step - 1;  x_fin(2) = this.last_x
		y_fin(:) = this.nf_y;  y_fin(3) = this.rcv_xy(2, 3, 2) + this.step - 1;  y_fin(1) = this.last_y

		do face = 1, 6
			do i = 1, 4
				! if(this.interp_factor(i) == 1) then
					do y = this.rcv_xy(2, i, 2), y_fin(i)
						do x = this.rcv_xy(2, i, 1), x_fin(i)

this.u_con(x, y, face) = metr.J_to_cube(1, 1, x, y, face) * this.lon_vel(x, y, face) + metr.J_to_cube(1, 2, x, y, face) * this.lat_vel(x, y, face)
this.v_con(x, y, face) = metr.J_to_cube(2, 2, x, y, face) * this.lat_vel(x, y, face) + metr.J_to_cube(2, 1, x, y, face) * this.lon_vel(x, y, face)

this.u_cov(x, y, face) = metr.G_tensor(1, 1, x, y) * this.u_con(x, y, face) + metr.G_tensor(1, 2, x, y) * this.v_con(x, y, face)
this.v_cov(x, y, face) = metr.G_tensor(2, 2, x, y) * this.v_con(x, y, face) + metr.G_tensor(2, 1, x, y) * this.u_con(x, y, face)

						end do
					end do
				! end if
			end do
		end do

	end subroutine




end module