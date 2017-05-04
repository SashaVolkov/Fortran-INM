module func_var

	use sphere_geometry, Only: geometry
	use interpolation, Only: interp
	use metrics, Only: metric

implicit none

	Private
	Public :: f_var

	Type f_var

		Real(8), Allocatable :: h_height(:, :, :)
		Real(8), Allocatable :: h_depth(:, :, :)
		Real(8), Allocatable :: starter(:, :, :, :)
		Real(8), Allocatable :: u_con(:, :, :)
		Real(8), Allocatable :: v_con(:, :, :)
		Real(8), Allocatable :: lon_vel(:, :, :)
		Real(8), Allocatable :: lat_vel(:, :, :)
		Real(8), Allocatable :: f(:, :, :)
		Real(8) height,  g, dt, delta_on_cube, u0
		Integer(4) step, dim, interp_factor(1:4), Neighbours_face(6, 4), grid_type, rescale, test
		Integer(4) ns_x, ns_y, nf_x, nf_y, first_x, first_y, last_x, last_y, snd_xy(6, 4, 2), rcv_xy(6, 4, 2)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Public :: equal => equal
		Procedure, Public :: start_conditions => start_conditions
		Procedure, Public :: interpolate => interpolate
		Procedure, Private :: Velocity_from_spherical_border => Velocity_from_spherical_border

	End Type


CONTAINS



	Subroutine init(this, metr, test)

		Class(f_var) :: this
		Class(metric) :: metr
		Integer(4), intent(in) :: test
		Integer(4) :: i, x, y, face

		this.g = 980616d-5;  this.test = test

		this.ns_x = metr.ns_xy(1);  this.ns_y = metr.ns_xy(2)
		this.nf_x = metr.nf_xy(1);  this.nf_y = metr.nf_xy(2)

		this.first_x = metr.first_x;  this.first_y = metr.first_y
		this.last_x = metr.last_x;  this.last_y = metr.last_y

		this.step = metr.step;  this.dim = metr.dim
		this.Neighbours_face = metr.Neighbours_face;  this.grid_type = metr.grid_type
		if ( test == 1 ) then
			this.height = 1000d0
		elseif ( test == 2 ) then
			this.height = 294d2/this.g;
		elseif ( test == 5 ) then
			this.height = 5960d0
		elseif ( test == 6 ) then
			this.height = 8d3
		end if

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
		Allocate(this.h_depth(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.starter(3, this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.u_con(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.v_con(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.lon_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.lat_vel(this.first_x:this.last_x, this.first_y:this.last_y, 6))
		Allocate(this.f(this.first_x:this.last_x, this.first_y:this.last_y, 6))

	end Subroutine



	Subroutine deinit(this)
		Class(f_var) :: this

		if (Allocated(this.h_height)) Deallocate(this.h_height)
		if (Allocated(this.h_depth)) Deallocate(this.h_depth)
		if (Allocated(this.starter)) Deallocate(this.starter)
		if (Allocated(this.u_con)) Deallocate(this.u_con)
		if (Allocated(this.v_con)) Deallocate(this.v_con)
		if (Allocated(this.lon_vel)) Deallocate(this.lon_vel)
		if (Allocated(this.lat_vel)) Deallocate(this.lat_vel)
		if (Allocated(this.f)) Deallocate(this.f)

	end Subroutine



	Subroutine equal(var_pr, var, metr)

		Class(f_var) :: var_pr, var
		Class(metric) :: metr
		Integer(4) :: face, x, y

		do face = 1, 6
			do y = var.ns_y, var.nf_y
				do x= var.ns_x, var.nf_x
		var_pr.h_height(x, y, face) = var.h_height(x, y, face)
		var_pr.u_con(x, y, face) = var.u_con(x, y, face)
		var_pr.v_con(x, y, face) = var.v_con(x, y, face)
		call metr.con_to_spherical(var_pr.u_con(x, y, face), var_pr.v_con(x, y, face), var_pr.lon_vel(x, y, face), var_pr.lat_vel(x, y, face), x, y, face)
				end do
			end do
		end do

	end Subroutine



	Subroutine start_conditions(this, metr, geom, omega_cor)

		Class(f_var) :: this
		Class(metric) :: metr
		Class(geometry) :: geom
		Real(8), intent(in) :: omega_cor
		Integer(4) dim, x, y, face, test
		Real(8) :: gh0, r, R_BIG, center(2), pi, u0, alpha, lon, lat, f, f_con_y, f_con_x, A_inv(2,2), G(2,2), h_s, h_s0
		Real(8) :: a, A_t, B_t, C_t, K, cos_lat, cos_lat2, p

		pi = 314159265358979323846d-20
		gh0 = 294d2

		dim = this.dim;  test = this.test
		this.lon_vel = 0d0;  this.u_con = 0d0;  this.v_con = 0d0;  this.lat_vel = 0d0

		center(:) = (/0d0, 3d0*pi/2d0/);  alpha = -pi/4d0
		this.u0 = 2d0*pi*metr.r_sphere/(12d0*24d0*60d0*60d0);  R_BIG = metr.r_sphere/3d0

		if ( test == 5 ) then
			center(:) = (/pi/6d0, 3d0*pi/2d0/);  alpha = 0d0
			this.u0 = 20d0;  R_BIG = pi/9d0;  h_s0 = 2000d0
		else if ( test == 6 ) then
			alpha = 0d0
		end if

		u0 = this.u0

		do face = 1, 6
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x

lat = metr.latlon_c(1,x,y,face)
lon = metr.latlon_c(2,x,y,face)

this.lon_vel(x, y, face) = u0*(dcos(lat)*dcos(alpha) + dcos(lon)*dsin(lat)*dsin(alpha))
this.lat_vel(x, y, face) = -u0*dsin(lon)*dsin(alpha)

if ( test == 1 ) then
	r = geom.dist(center(:),metr.latlon_c(1:2,x,y,face))  ! cosine bell
	this.h_height(x, y, face) = (this.height/2d0)*(1d0 + dcos(pi*r/R_BIG))
	if (r >= R_BIG) this.h_height(x, y, face) = 0d0
else if ( test == 2 ) then
	this.h_height(x, y, face) = this.height - (metr.r_sphere*omega_cor*u0 + 5d-1*u0*u0)*((dsin(lat)*dcos(alpha) - dcos(lat)*dcos(lon)*dsin(alpha))**2)/this.g
else if ( test == 5 ) then
	r = sqrt(min( R_BIG**2, (center(2) - lon)**2 + (center(1) + lat)**2))
	h_s = h_s0*(1d0 - r/R_BIG)
	this.h_depth(x, y, face) = h_s0 - h_s
	this.h_height(x, y, face) = this.height - h_s0
else if ( test == 6 ) then
	a = metr.r_sphere;  K = 789d-9;  cos_lat = dcos(lat);  cos_lat2 = cos_lat**2;  R_BIG = 4;  p = (R_BIG + 1d0)*cos_lat2

	A_t = (K/2d0)*(2d0*omega_cor + K)*cos_lat2 + 	25d-2 * (K**2) * (dcos(lat*(p + (2*R_BIG**2 - R_BIG - 2d0) - (2*R_BIG**2)/cos_lat2)))**(2d0*R_BIG)
	B_t = (2d0*(omega_cor + K)*K/((R_BIG+1)*(R_BIG+2)))*(cos_lat**R_BIG)*((R_BIG**2 + 2*R_BIG + 2) - ((R_BIG+1)**2)*cos_lat2)
	C_t = 25d-2*(K**2)*(cos_lat2**R_BIG)*((R_BIG+1)*cos_lat2 - (R_BIG+2))

	this.h_height(x, y, face) = this.height + (a**2)*(A_t + B_t*dcos(R_BIG*lon) + C_t*dcos(2d0*R_BIG*lon))/this.g
	this.lon_vel(x, y, face) = a*K*cos_lat + a*K*(cos_lat**(R_BIG-1))*(R_BIG*(1-cos_lat2) - cos_lat2)*dcos(R_BIG*lon)
	this.lat_vel(x, y, face) = -a*K*R_BIG*(cos_lat**(R_BIG-1))*dsin(lat)*dsin(R_BIG*lon)
end if



this.f(x, y, face) = (2d0*omega_cor)*(dsin(lat)*dcos(alpha) - dcos(lat)*dcos(lon)*dsin(alpha))
call metr.spherical_to_con(this.lon_vel(x, y, face), this.lat_vel(x, y, face), this.u_con(x, y, face), this.v_con(x, y, face), x, y, face)
				end do
			end do
		end do

		this.starter(1,:,:,:) = this.h_height
		this.starter(2,:,:,:) = this.lon_vel
		this.starter(3,:,:,:) = this.lat_vel

	end Subroutine




	Subroutine interpolate(this, i, metr)
		Class(f_var) :: this
		Class(interp) :: i
		Class(metric) :: metr

		call i.Lagrange(this.h_height, this.interp_factor)
		call i.Lagrange(this.lat_vel, this.interp_factor)
		call i.Lagrange(this.lon_vel, this.interp_factor)
		call this.Velocity_from_spherical_border(metr)

	end Subroutine



	Subroutine Velocity_from_spherical_border(this, metr)
		Class(f_var) :: this
		Class(metric) :: metr
		Real(8) :: G(2,2), A_inv(2,2)
		Integer(4) :: x, y, face, i, x_fin(4), y_fin(4), x_st(4), y_st(4)

		x_st(1) = this.first_x;  y_st(1) = this.nf_y + 1;  x_fin(1) = this.last_x ;  y_fin(1) = this.last_y
		x_st(2) = this.nf_x+1;  y_st(2) = this.ns_y;  x_fin(2) = this.last_x ;  y_fin(2) = this.nf_y
		x_st(3) = this.first_x;  y_st(3) = this.first_y;  x_fin(3) = this.last_x ;  y_fin(3) = this.ns_y-1
		x_st(4) = this.first_x;  y_st(4) = this.ns_y;  x_fin(4) = this.ns_x-1 ;  y_fin(4) = this.nf_y


		do face = 1, 6
			do i = 1, 4
				do y = y_st(i), y_fin(i)
					do x = x_st(i), x_fin(i)

call metr.spherical_to_con(this.lon_vel(x, y, face), this.lat_vel(x, y, face), this.u_con(x, y, face), this.v_con(x, y, face), x, y, face)

					end do
				end do
			end do
		end do

	end Subroutine




end module