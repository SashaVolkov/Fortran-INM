module metrics

	use grid_var, Only: g_var

implicit none

	Private
	Public :: metric

	Type metric

		Real(8), Allocatable :: G_sqr(:, :)
		Real(8), Allocatable :: Christoffel_x1(:, :, :, :)
		Real(8), Allocatable :: Christoffel_x2(:, :, :, :)
		Real(8), Allocatable :: G_tensor(:, :, :, :)
		Real(8), Allocatable :: G_inverse(:, :, :, :)
		Real(8), Allocatable :: rho(:, :)
		Real(8), Allocatable :: Tr_to_sph(:, :, :, :, :)
		Real(8), Allocatable :: Tr_to_cube(:, :, :, :, :)
		Real(8), Allocatable :: latlon_c(:, :, :, :)
		Real(8), Allocatable :: cube_coord_c(:, :, :)
		Real(8), Allocatable :: S_cor(:, :, :, :, :)

		Real(8) :: r_sphere, delta_on_cube
		Integer(4) dim, step, rescale, ns_xy(2), nf_xy(2), snd_xy(6, 4, 2), rcv_xy(6, 4, 2)
		Integer(4) first_x, first_y, last_x, last_y, grid_type, Neighbours_face(6, 4)

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: deinit => deinit
		Procedure, Public :: define => define
		Procedure, Private :: alloc => alloc

		Procedure, Private :: transf_matrix_equiang => transf_matrix_equiang
		Procedure, Private :: metric_tensor_equiang => metric_tensor_equiang

		Procedure, Private :: transf_matrix_conf => transf_matrix_conf
		Procedure, Private :: metric_tensor_conf => metric_tensor_conf

		Procedure, Private :: hem_of_face => hem_of_face
		Procedure, Private :: partial => partial

		Procedure, Public :: cov_to_con => cov_to_con
		Procedure, Public :: con_to_cov => con_to_cov

		Procedure, Public :: spherical_to_con => spherical_to_con
		Procedure, Public :: con_to_spherical => con_to_spherical

	End Type


CONTAINS


	Subroutine init(this, grid)

		Class(metric) :: this
		Class(g_var) :: grid
		Integer(4) :: face, x, y


		this.dim = grid.dim;  this.step = grid.step
		this.delta_on_cube = grid.delta_on_cube

		this.ns_xy(:) = grid.ns_xy(:);  this.nf_xy(:) = grid.nf_xy(:);
		this.first_x = grid.first_x;  this.first_y = grid.first_y
		this.last_x = grid.last_x;  this.last_y = grid.last_y

		call this.alloc()

		this.r_sphere = grid.r_sphere;  this.grid_type = grid.grid_type
		this.snd_xy = grid.snd_xy;  this.rcv_xy = grid.rcv_xy
		this.Neighbours_face = grid.Neighbours_face;  this.rescale = grid.rescale

		this.cube_coord_c = grid.cube_coord_c
		this.latlon_c = grid.latlon_c
		call this.define()

	end Subroutine


	Subroutine define(this)

		Class(metric) :: this

		if(this.grid_type == 0) then ! 0 - conformal, 1 - equiangular
			call this.metric_tensor_conf()
		else if(this.grid_type == 1)then
			call this.metric_tensor_equiang()
		end if


	end Subroutine



	Subroutine alloc(this)
		Class(metric) :: this
		Integer(4) f_x, f_y, l_x, l_y, dim, step, f, l

		f_x = this.first_x;  l_x = this.last_x;  f_y = this.first_y;  l_y = this.last_y
		dim = this.dim;  step = this.step;  f = 1-2*step; l = 2*dim + 2*step

		Allocate(this.G_sqr(f_x:l_x , f_y:l_y))
		Allocate(this.G_tensor(2, 2, f_x:l_x , f_y:l_y))
		Allocate(this.Christoffel_x1(3, 3, f_x:l_x , f_y:l_y))
		Allocate(this.Christoffel_x2(3, 3, f_x:l_x , f_y:l_y))
		Allocate(this.G_inverse(2, 2, f_x:l_x , f_y:l_y))
		Allocate(this.rho(f_x:l_x , f_y:l_y))
		Allocate(this.Tr_to_sph(2, 2, f_x:l_x , f_y:l_y, 6))
		Allocate(this.Tr_to_cube(2, 2, f_x:l_x , f_y:l_y, 6))
		Allocate(this.cube_coord_c(2, f:l , f:l))
		Allocate(this.latlon_c(2, f:l , f:l, 6))
		Allocate(this.S_cor(2, 2, f:l , f:l, 6))

	end Subroutine



	Subroutine deinit(this)
		Class(metric) :: this

		if (Allocated(this.G_sqr)) Deallocate(this.G_sqr)
		if (Allocated(this.G_tensor)) Deallocate(this.G_tensor)
		if (Allocated(this.Christoffel_x1)) Deallocate(this.Christoffel_x1)
		if (Allocated(this.Christoffel_x2)) Deallocate(this.Christoffel_x2)
		if (Allocated(this.G_inverse)) Deallocate(this.G_inverse)
		if (Allocated(this.rho)) Deallocate(this.rho)
		if (Allocated(this.Tr_to_sph)) Deallocate(this.Tr_to_sph)
		if (Allocated(this.Tr_to_cube)) Deallocate(this.Tr_to_cube)
		if (Allocated(this.latlon_c)) Deallocate(this.latlon_c)
		if (Allocated(this.cube_coord_c)) Deallocate(this.cube_coord_c)

	end Subroutine



	Subroutine metric_tensor_equiang(this)   ! Ullrich phd thesis Appendices
		Class(metric) :: this
		Real(8) x_1, x_2, g_coef, g_inv_coef, delta, rho
		Integer(4) x, y, face

		this.Christoffel_x1 = 0d0
		this.Christoffel_x2 = 0d0
		call this.transf_matrix_equiang()

		do y = this.first_y, this.last_y
			do x = this.first_x, this.last_x

			x_1 = dtan(this.cube_coord_c(1, x, y))
			x_2 = dtan(this.cube_coord_c(2, x, y))

			rho = 1d0 + x_1**2 + x_2**2
			delta = dsqrt(rho)
			this.G_sqr(x, y) = (1d0 + x_1**2)*(1d0 + x_2**2)/(delta*rho)
			g_coef = (1d0 + x_1**2)*(1d0 + x_2**2)/(rho**2)
			g_inv_coef = (rho)/((1d0 + x_1**2)*(1d0 + x_2**2))

			this.G_tensor(1, 1, x, y) = g_coef * (1d0 + x_1**2)
			this.G_tensor(1, 2, x, y) = g_coef * (- x_1*x_2)
			this.G_tensor(2, 1, x, y) = this.G_tensor(1, 2, x, y)
			this.G_tensor(2, 2, x, y) = g_coef * (1d0 + x_2**2)

			this.G_inverse(1, 1, x, y) = g_inv_coef * (1d0 + x_2**2)
			this.G_inverse(1, 2, x, y) = g_inv_coef * (x_1*x_2)
			this.G_inverse(2, 1, x, y) = this.G_inverse(1, 2, x, y)
			this.G_inverse(2, 2, x, y) = g_inv_coef * (1d0 + x_1**2)

			rho = rho*this.r_sphere
			this.Christoffel_x1(1, 1, x, y) = 2d0*x_1*(x_2**2)/rho   ! Ullrich PhD thesis
			this.Christoffel_x1(1, 2, x, y) = -x_2*(1 + x_2**2)/rho
			this.Christoffel_x1(2, 1, x, y) = -x_2*(1 + x_2**2)/rho

			this.Christoffel_x2(2, 2, x, y) = 2d0*x_2*(x_1**2)/rho
			this.Christoffel_x2(1, 2, x, y) = -x_1*(1 + x_1**2)/rho
			this.Christoffel_x2(2, 1, x, y) = -x_1*(1 + x_1**2)/rho

			end do
		end do

	end Subroutine



	Subroutine transf_matrix_equiang(this)   ! Ullrich phd thesis Appendix G.4
		Class(metric) :: this
		Real(8) :: x_1, x_2, g_coef, s(6), cos_theta, delta, omega_cor, f, lat, A(2,2), det, pi
		Integer(4) :: x, y, face
		s(1) = - 1d0;  s(6) = 1d0;  omega_cor = 7292d-8
		this.Tr_to_cube = 0d0;  this.Tr_to_sph = 0d0;  pi = 314159265358979323846d-20


		do face = 2,5
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					lat = this.latlon_c(1, x, y, face)
					x_1 = dtan(this.cube_coord_c(1, x, y))
					x_2 = dtan(this.cube_coord_c(2, x, y))
					delta = 1d0 + x_1**2 + x_2**2

					this.Tr_to_sph(1,1,x,y,face) = dcos(lat)
					this.Tr_to_sph(1,2,x,y,face) = 0d0
					this.Tr_to_sph(2,1,x,y,face) = - x_1*x_2*dsqrt(1d0 + x_1**2)/(delta)
					this.Tr_to_sph(2,2,x,y,face) = ((1d0 + x_2**2)*dsqrt(1d0 + x_1**2))/(delta)

				end do
			end do
		end do


		do face = 1, 6, 5
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					lat = this.latlon_c(1, x, y, face)
					x_1 = dtan(this.cube_coord_c(1, x, y))
					x_2 = dtan(this.cube_coord_c(2, x, y))
					delta = 1d0 + x_1**2 + x_2**2

					this.Tr_to_sph(1,1,x,y,face) = -s(face)*x_2*(1d0 + x_1**2)*dcos(lat)/(x_2**2 + x_1**2)
					this.Tr_to_sph(1,2,x,y,face) = s(face)*x_1*(1d0 + x_2**2)*dcos(lat)/(x_2**2 + x_1**2)
					this.Tr_to_sph(2,1,x,y,face) = - s(face)*x_1*(1d0 + x_1**2)/((delta)*dsqrt(x_2**2 + x_1**2))
					this.Tr_to_sph(2,2,x,y,face) = - s(face)*x_2*(1d0 + x_2**2)/((delta)*dsqrt(x_2**2 + x_1**2))

				end do
			end do
		end do


		do face = 1, 6
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x

				A = this.Tr_to_sph(:,:,x,y,face)
				det = A(1,1)*A(2,2) - A(2,1)*A(1,2)

				this.Tr_to_cube(1,1,x,y,face) = A(2,2)/det
				this.Tr_to_cube(1,2,x,y,face) = -A(1,2)/det
				this.Tr_to_cube(2,1,x,y,face) = -A(2,1)/det
				this.Tr_to_cube(2,2,x,y,face) = A(1,1)/det

				end do
			end do
		end do

	end Subroutine





	Subroutine metric_tensor_conf(this) ! Tranformation matrix /= Jacobi matrix ! Tranformation matrix == J^T
		Class(metric) :: this
		Integer(4) x, y, i, k, dim
		Real(8) :: A(2,2), G(2,2), det, cos_theta

		dim = this.dim
		call this.transf_matrix_conf()

		do y = this.first_y, this.last_y
			do x = this.first_x, this.last_x

			A(:,:) = this.Tr_to_sph(:, :, x, y, 2)
			cos_theta = dcos(this.latlon_c(1, x, y, 2))

			this.G_tensor(1, 1, x, y) = (A(1, 1)**2)*(cos_theta**2) + (A(2, 1)**2)
			this.G_tensor(1, 2, x, y) = (A(1, 1)*A(1,2))*(cos_theta**2) + A(2, 1)*A(2,2)
			this.G_tensor(2, 1, x, y) = this.G_tensor(1, 2, x, y)
			this.G_tensor(2, 2, x, y) = (A(1, 2)**2)*(cos_theta**2) + A(2, 2)**2

			G(:,:) = this.G_tensor(:, :, x, y)
			det = G(1,1)*G(2,2) - G(2,1)*G(1,2)
			this.G_sqr(x, y) = dsqrt(det)
			if(this.G_sqr(x, y) == 0d0) print*, "Divide by zero:", x, y

			this.G_inverse(1, 1, x, y) = this.G_tensor(2, 2, x, y)/det
			this.G_inverse(1, 2, x, y) = -this.G_tensor(1, 2, x, y)/det
			this.G_inverse(2, 1, x, y) = -this.G_tensor(2, 1, x, y)/det
			this.G_inverse(2, 2, x, y) = this.G_tensor(1, 1, x, y)/det

			end do
		end do


	end Subroutine



	Subroutine transf_matrix_conf(this)
		Class(metric) :: this
		Real(8) x_1, x_2, g_coef, delta, temp(-this.step:this.step), A(2,2), det
		Integer(4) x, y, dim, i, k, face, step

		dim = this.dim;  step = this.step
		delta = (1d0/dble(dim))

		call this.hem_of_face(this.latlon_c(:, :, :, 1:6), 1)

		do face = 1, 6
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					temp = this.latlon_c(2, x-step:x+step, y, face)
					this.Tr_to_sph(1,1,x,y,face) = this.partial(temp, delta)

					temp = this.latlon_c(1, x-step:x+step, y, face)
					this.Tr_to_sph(2,1,x,y,face) = this.partial(temp, delta)
				end do
			end do
		end do

		call this.hem_of_face(this.latlon_c(:, :, :, 1:6), 2)

		do face = 1, 6
			do x = this.first_x, this.last_x
				do y = this.first_y, this.last_y
					temp = this.latlon_c(2, x, y-step:y+step, face)
					this.Tr_to_sph(1,2,x,y,face) = this.partial(temp, delta)

					temp = this.latlon_c(1, x, y-step:y+step, face)
					this.Tr_to_sph(2,2,x,y,face) = this.partial(temp, delta)

					A = this.Tr_to_sph(:,:,x,y,face)
					det = A(1,1)*A(2,2) - A(2,1)*A(1,2)

					this.Tr_to_cube(1,1,x,y,face) = A(2,2)/det
					this.Tr_to_cube(1,2,x,y,face) = -A(1,2)/det
					this.Tr_to_cube(2,1,x,y,face) = -A(2,1)/det
					this.Tr_to_cube(2,2,x,y,face) = A(1,1)/det

				end do
			end do
		end do



	end Subroutine



	Subroutine hem_of_face(this, cubic, x_or_y)
		Class(metric) :: this
		Real(8), intent(inout) :: cubic(1:2, 1 - 2*this.step:2*this.dim + 2*this.step, 1 - 2*this.step:2*this.dim + 2*this.step, 6)
		Integer(4), intent(in) :: x_or_y
		Integer(4) dim, f, l, x, y, face, i, j, k, step

		dim = this.dim;  step = this.step
		f = 1 - 2*step; l = 2*dim + 2*step

		do i = 1, 2*step
			do j = 1, 2*dim

				cubic(:,j,2*dim+i,2) = cubic(:,j,i,6);      cubic(:,j,1-i,6) = cubic(:,j,2*dim+1-i,2)
				cubic(:,2*dim+i,j,2) = cubic(:,i,j,3);      cubic(:,1-i,j,3) = cubic(:,2*dim+1-i,j,2)
				cubic(:,j,1-i,2) = cubic(:,j,2*dim+1-i,1);  cubic(:,j,2*dim+i,1) = cubic(:,j,i,2)
				cubic(:,1-i,j,2) = cubic(:,2*dim+1-i,j,5);  cubic(:,2*dim+i,j,5) = cubic(:,i,j,2)

			end do
		end do

		do i = 1, 2*step
			do j = 1, 2*dim
			k = 2*dim + 1 -j

				cubic(:,j,2*dim+i,4) = cubic(:,k,2*dim+1-i,6);     cubic(:,k,2*dim+i,6) = cubic(:,j,2*dim+1-i,4)
				cubic(:,2*dim+i,j,4) = cubic(:,i,j,5);             cubic(:,1-i,j,5) = cubic(:,2*dim+1-i,j,4)
				cubic(:,j,1-i,4) = cubic(:,k,i,1);                 cubic(:,k,1-i,1) = cubic(:,j,i,4)
				cubic(:,1-i,j,4) = cubic(:,2*dim+1-i,j,3);         cubic(:,2*dim+i,j,3) = cubic(:,i,j,4)

			end do
		end do

		do i = 1, 2*step
			do j = 1, 2*dim
			k = 2*dim + 1 -j

				cubic(:,j,2*dim+i,3) = cubic(:,2*dim+1-i,j,6);     cubic(:,2*dim+i,j,6) = cubic(:,j,2*dim+1-i,3)
				cubic(:,j,1-i,3) = cubic(:,2*dim+1-i,k,1);         cubic(:,2*dim+i,k,1) = cubic(:,j,i,3)

			end do
		end do

		do i = 1, 2*step
			do j = 1, 2*dim
			k = 2*dim + 1 -j

				cubic(:,j,2*dim+i,5) = cubic(:,i,k,6);     cubic(:,1-i,k,6) = cubic(:,j,2*dim+1-i,5)
				cubic(:,j,1-i,5) = cubic(:,i,j,1);         cubic(:,1-i,j,1) = cubic(:,j,i,5)

			end do
		end do

		do face = 1, 6
		do i = 1, 2*step
			do j = 1, 2*step
				if( x_or_y == 1) then
					cubic(:,1-i, 1-j, face) = cubic(:,1-j, i, face)
					cubic(:,1-i, 2*dim + j, face) = cubic(:,1-j, 2*dim+1-i, face)

					cubic(:,2*dim + i, 1-j, face) = cubic(:,2*dim+j, i, face)
					cubic(:,2*dim + i, 2*dim + j, face) = cubic(:,2*dim+j, 2*dim+1-i, face)
				else if(x_or_y == 2) then
					cubic(:,1-i, 1-j, face) = cubic(:,j, 1-i, face)
					cubic(:,1-i, 2*dim + j, face) = cubic(:,j, 2*dim+i, face)

					cubic(:,2*dim + i, 1-j, face) = cubic(:,2*dim+1-j, 1-i, face)
					cubic(:,2*dim + i, 2*dim + j, face) = cubic(:,2*dim+1-j, 2*dim+i, face)
				end if
			end do
		end do
		end do

	end Subroutine



	Subroutine cov_to_con(this, u_cov, v_cov, u_con, v_con, x, y)
		Class(metric) :: this
		Integer(4), intent(in) :: x, y
		Real(8), intent(in) :: u_cov, v_cov
		Real(8), intent(out) :: u_con, v_con
		Real(8) :: G_inv(2,2)

G_inv(:,:) = this.G_inverse(:,:, x, y)
u_con = G_inv(1, 1) * u_cov + G_inv(1, 2) * v_cov
v_con = G_inv(2, 2) * v_cov + G_inv(2, 1) * u_cov

	end Subroutine



	Subroutine con_to_cov(this, u_con, v_con, u_cov, v_cov, x, y)
		Class(metric) :: this
		Integer(4), intent(in) :: x, y
		Real(8), intent(in) :: u_con, v_con
		Real(8), intent(out) :: u_cov, v_cov
		Real(8) :: G(2,2)

G(:,:) = this.G_tensor(:,:, x, y)
u_cov = G(1, 1) * u_con + G(1, 2) * v_con
v_cov = G(2, 2) * v_con + G(2, 1) * u_con

	end Subroutine


		Subroutine con_to_spherical(this, u_con, v_con, lon_vel, lat_vel, x, y, face)
			Class(metric) :: this
			Integer(4), intent(in) :: x, y, face
			Real(8), intent(in) :: u_con, v_con
			Real(8), intent(out) :: lon_vel, lat_vel
			Real(8) :: vel_x_contr, vel_y_contr, A(2,2)

	A(:,:) = this.Tr_to_sph(:,:, x, y, face)
	lon_vel = A(1, 1) * u_con + A(1, 2) * v_con
	lat_vel = A(2, 2) * v_con + A(2, 1) * u_con

		end Subroutine



		Subroutine spherical_to_con(this, lon_vel, lat_vel, u_con, v_con, x, y, face)
			Class(metric) :: this
			Integer(4), intent(in) :: x, y, face
			Real(8), intent(in) :: lon_vel, lat_vel
			Real(8), intent(out) :: u_con, v_con
			Real(8) :: A_inv(2,2)

	A_inv(:,:) = this.Tr_to_cube(:,:, x, y, face)
	u_con = A_inv(1, 1) * lon_vel + A_inv(1, 2) * lat_vel
	v_con = A_inv(2, 2) * lat_vel + A_inv(2, 1) * lon_vel

		end Subroutine



	Real(8) function partial(this, fun, h)
		Class(metric) :: this
		Real(8), intent(in) :: fun(-this.step:this.step), h
		Real(8) A , B, C, D, E, F, G, I

		if(this.step == 2) then
			A = 2d0/(3d0*h);  B = - 2d0/(3d0*h);  C = - 1d0/(12d0*h);  D = 1d0/(12d0*h)
			partial = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2)
		else if(this.step == 3) then
			A = 3d0/(4d0*h);  B = - 3d0/(4d0*h);  C = - 3d0/(20d0*h);  D = 3d0/(20d0*h)
			E = 1d0/(60d0*h);  F = - 1d0/(60d0*h)
			partial = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) + E*fun(3) + F*fun(-3)
		else if(this.step == 4) then
			A = 4d0/(5d0*h);  B = - 4d0/(5d0*h);  C = - 1d0/(5d0*h);  D = 1d0/(5d0*h)
			E = 4d0/(105d0*h);  F = - 4d0/(105d0*h); G = - 1d0/(280d0*h); I = 1d0/(280d0*h)

			partial = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) + E*fun(3) + F*fun(-3) + G*fun(4) + I*fun(-4)
		end if

	end function



end module