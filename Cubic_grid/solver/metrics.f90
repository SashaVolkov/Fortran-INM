module metrics

	use parallel_cubic, Only: parallel

implicit none

	Private
	Public :: metric

	Type metric

		Real(8), Allocatable :: G_sqr(:, :)
		Real(8), Allocatable :: G_tensor(:, :, :, :)
		Real(8), Allocatable :: G_inverse(:, :, :, :)
		Real(8), Allocatable :: rho(:, :)
		Real(8), Allocatable :: J_to_sph(:, :, :, :, :)
		Real(8), Allocatable :: J_to_cube(:, :, :, :, :)
		Real(8), Allocatable :: latlon_c(:, :, :, :)
		Real(8), Allocatable :: cube_coord_c(:, :, :)

		real(8) :: r_sphere
		integer(4) dim, step, rescale, ns_xy(2), nf_xy(2)
		integer(4) first_x, first_y, last_x, last_y, grid_type

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
		Procedure, Private :: partial_c4 => partial_c4


	End Type


CONTAINS


	subroutine init(this, paral)

		Class(metric) :: this
		Class(parallel) :: paral


		this.dim = paral.dim;  this.step = paral.step

		this.ns_xy(:) = paral.ns_xy(:);  this.nf_xy(:) = paral.nf_xy(:);
		this.first_x = paral.first_x;  this.first_y = paral.first_y
		this.last_x = paral.last_x;  this.last_y = paral.last_y

		call this.alloc()

	end subroutine


	subroutine define(this, grid_type)

		Class(metric) :: this
		integer(4), intent(in) :: grid_type

		this.grid_type = grid_type

		if(grid_type == 0) then ! 0 - conformal, 1 - equiangular
			call this.metric_tensor_conf()
		else if(grid_type == 1)then
			call this.metric_tensor_equiang()
		end if


	end subroutine



	subroutine alloc(this)
		Class(metric) :: this
		integer(4) f_x, f_y, l_x, l_y, dim, step, f, l

		f_x = this.first_x;  l_x = this.last_x;  f_y = this.first_y;  l_y = this.last_y
		dim = this.dim;  step = this.step;  f = 1-step; l = 2*dim + step

		Allocate(this.G_sqr(f:l , f:l))
		Allocate(this.G_tensor(2, 2, f:l , f:l))
		Allocate(this.G_inverse(2, 2, f:l , f:l))
		Allocate(this.rho(f:l , f:l))
		Allocate(this.J_to_sph(2, 2, f:l , f:l, 6))
		Allocate(this.J_to_cube(2, 2, f:l , f:l, 6))
		Allocate(this.cube_coord_c(1:2, f:l , f:l))
		Allocate(this.latlon_c(1:2, f:l , f:l, 1:6))

	end subroutine



	subroutine deinit(this)
		Class(metric) :: this

		if (Allocated(this.G_sqr)) Deallocate(this.G_sqr)
		if (Allocated(this.G_tensor)) Deallocate(this.G_tensor)
		if (Allocated(this.G_inverse)) Deallocate(this.G_inverse)
		if (Allocated(this.rho)) Deallocate(this.rho)
		if (Allocated(this.J_to_sph)) Deallocate(this.J_to_sph)
		if (Allocated(this.J_to_cube)) Deallocate(this.J_to_cube)
		if (Allocated(this.latlon_c)) Deallocate(this.latlon_c)
		if (Allocated(this.cube_coord_c)) Deallocate(this.cube_coord_c)

	end subroutine



	subroutine metric_tensor_equiang(this)   ! Ullrich phd thesis Appendices
		Class(metric) :: this
		real(8) x_1, x_2, g_coef, g_inv_coef
		integer(4) x, y, face

		call this.transf_matrix_equiang()

		do y = this.first_y, this.last_y
			do x = this.first_x, this.last_x

			x_1 = dtan(this.cube_coord_c(1, x, y))
			x_2 = dtan(this.cube_coord_c(2, x, y))

			this.rho(x, y) = dsqrt(1d0 + x_1**2 + x_2**2)
			this.G_sqr(x, y) = (1d0 + x_1**2)*(1d0 + x_2**2)/(this.rho(x, y)**3)
			g_coef = (1d0 + x_1**2)*(1d0 + x_2**2)/(this.rho(x, y)**4)
			g_inv_coef = (this.rho(x, y)**2)/((1d0 + x_1**2)*(1d0 + x_2**2))

			this.G_tensor(1, 1, x, y) = g_coef * (1d0 + x_1**2)
			this.G_tensor(1, 2, x, y) = g_coef * (- x_1*x_2)
			this.G_tensor(2, 1, x, y) = this.G_tensor(1, 2, x, y)
			this.G_tensor(2, 2, x, y) = g_coef * (1d0 + x_2**2)

			this.G_inverse(1, 1, x, y) = g_inv_coef * (1d0 + x_2**2)
			this.G_inverse(1, 2, x, y) = g_inv_coef * (x_1*x_2)
			this.G_inverse(2, 1, x, y) = this.G_inverse(1, 2, x, y)
			this.G_inverse(2, 2, x, y) = g_inv_coef * (1d0 + x_1**2)

			end do
		end do

	end subroutine



	subroutine transf_matrix_equiang(this)   ! Ullrich phd thesis Appendix G.4
		Class(metric) :: this
		real(8) x_1, x_2, g_coef, s(6), cos_theta, delta
		integer(4) x, y, face
		s(1) = - 1d0;  s(6) = 1d0


		do face = 2,5
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					x_1 = dtan(this.cube_coord_c(1, x, y))
					x_2 = dtan(this.cube_coord_c(2, x, y))
					delta = dsqrt(1d0 + x_1**2 + x_2**2)
					cos_theta = dcos(this.latlon_c(1, x, y, face))

					this.J_to_cube(1,1,x,y,face) = 1d0
					this.J_to_cube(1,2,x,y,face) = 0d0
					this.J_to_cube(2,1,x,y,face) = x_1*x_2/(1 + x_2**2)
					this.J_to_cube(2,2,x,y,face) = (delta**2)/((1d0 + x_2**2)*dsqrt(1d0 + x_1**2))

					this.J_to_sph(1,1,x,y,face) = 1d0
					this.J_to_sph(1,2,x,y,face) = 0d0
					this.J_to_sph(2,1,x,y,face) = - x_1*x_2*dsqrt(1d0 + x_1**2)/(delta**2)
					this.J_to_sph(2,2,x,y,face) = ((1d0 + x_2**2)*dsqrt(1d0 + x_1**2))/(delta**2)
				end do
			end do
		end do


		do face = 1, 6, 5
			do y = this.first_y, this.last_y
				do x = this.first_x, this.last_x
					x_1 = dtan(this.cube_coord_c(1, x, y))
					x_2 = dtan(this.cube_coord_c(2, x, y))
					delta = dsqrt(1d0 + x_1**2 + x_2**2)
					cos_theta = dcos(this.latlon_c(1, x, y, face))


					this.J_to_cube(1,1,x,y,face) = -s(face)*x_2/(1d0 + x_1**2)
					this.J_to_cube(1,2,x,y,face) = -s(face)*(delta**2)*x_1/((1d0 + x_1**2)*dsqrt(x_2**2 + x_1**2))
					this.J_to_cube(2,1,x,y,face) = s(face)*x_1/(1d0 + x_2**2)
					this.J_to_cube(2,2,x,y,face) = -s(face)*(delta**2)*x_2/((1d0 + x_2**2)*dsqrt(x_2**2 + x_1**2))

					this.J_to_sph(1,1,x,y,face) = -s(face)*x_2*(1d0 + x_1**2)/(x_2**2 + x_1**2)
					this.J_to_sph(1,2,x,y,face) = s(face)*x_1*(1d0 + x_2**2)/(x_2**2 + x_1**2)
					this.J_to_sph(2,1,x,y,face) = - s(face)*x_1*(1d0 + x_1**2)/((delta**2)*dsqrt(x_2**2 + x_1**2))
					this.J_to_sph(2,2,x,y,face) = - s(face)*x_2*(1d0 + x_2**2)/((delta**2)*dsqrt(x_2**2 + x_1**2))

				end do
			end do
		end do

	end subroutine



	subroutine metric_tensor_conf(this)
		Class(metric) :: this
		integer(4) x, y, i, k, dim
		real(8) :: J(2,2), cos_theta

		dim = this.dim
		call this.transf_matrix_conf()

		do y = 1-this.step, 2*dim+this.step
			do x = 1-this.step, 2*dim+this.step
			J(:,:) = this.J_to_sph(:, :, x, y, 2)
			cos_theta = dcos(this.latlon_c(1, x, y, 2))

			this.G_tensor(1, 1, x, y) = (J(1, 1)**2)*(cos_theta**2) + J(2, 1)**2
			this.G_tensor(1, 2, x, y) = (J(1, 1)*J(1,2))*(cos_theta**2) + J(2, 1)*J(2,2)
			this.G_tensor(2, 1, x, y) = this.G_tensor(1, 2, x, y)
			this.G_tensor(2, 2, x, y) = (J(1, 2)**2)*(cos_theta**2) + J(2, 2)**2

			this.G_inverse(1, 1, x, y) = 1d0/this.G_tensor(1, 1, x, y)
			this.G_inverse(1, 2, x, y) = - this.G_tensor(1, 2, x, y)
			this.G_inverse(2, 1, x, y) = - this.G_tensor(2, 1, x, y)
			this.G_inverse(2, 2, x, y) = 1d0/this.G_tensor(2, 2, x, y)

			! print *, real(this.G_tensor(:,:,x,y),4), x, y
			end do
		end do

		do i = 1, this.step
			do k = 1, 2*dim
				this.G_tensor(:, :, k, 2*dim + i) = this.G_tensor(:, :, k, i)
				this.G_tensor(:, :, 2*dim + i, k) = this.G_tensor(:, :, i, k)
				this.G_tensor(:, :, k, 1 - i) = this.G_tensor(:, :, k, 2*dim - i + 1)
				this.G_tensor(:, :, 1 - i, k) = this.G_tensor(:, :, 2*dim - i + 1, k)

				this.G_inverse(:, :, k, 2*dim + i) = this.G_inverse(:, :, k, i)
				this.G_inverse(:, :, 2*dim + i, k) = this.G_inverse(:, :, i, k)
				this.G_inverse(:, :, k, 1 - i) = this.G_inverse(:, :, k, 2*dim - i + 1)
				this.G_inverse(:, :, 1 - i, k) = this.G_inverse(:, :, 2*dim - i + 1, k)
			end do
		end do

		do y = 1 - this.step, 2*dim + this.step
			do x = 1 - this.step, 2*dim + this.step
				this.G_sqr(x, y) = dsqrt(this.G_tensor(1, 1, x, y) * this.G_tensor(2, 2, x, y) - this.G_tensor(2, 1, x, y) * this.G_tensor(1, 2, x, y))
			end do
		end do


	end subroutine



	subroutine transf_matrix_conf(this)
		Class(metric) :: this
		real(8) x_1, x_2, g_coef, delta, temp(-2:2)
		integer(4) x, y, dim, i, k, face

		dim = this.dim
		delta = (this.cube_coord_c(1, dim, dim) - this.cube_coord_c(1, dim-1, dim))

		call this.hem_of_face(this.latlon_c(1, :, :, 1:6))
		call this.hem_of_face(this.latlon_c(2, :, :, 1:6))

		do face = 1, 6
			do y = 1, 2*dim
				do x = 1, 2*dim

					temp = this.latlon_c(2, x-2:x+2, y, face)
					this.J_to_sph(1,1,x,y,face) = this.partial_c4(temp, delta)
					this.J_to_cube(1,1,x,y,face) = 2d0*delta/(temp(1) - temp(-1))

					temp = this.latlon_c(2, x, y-2:y+2, face)
					this.J_to_sph(1,2,x,y,face) = this.partial_c4(temp, delta)
					this.J_to_cube(2,1,x,y,face) = 2d0*delta/(temp(1) - temp(-1))

					temp = this.latlon_c(1, x-2:x+2, y, face)
					this.J_to_sph(2,1,x,y,face) = this.partial_c4(temp, delta)
					this.J_to_cube(1,2,x,y,face) = 2d0*delta/(temp(1) - temp(-1))

					temp = this.latlon_c(1, x, y-2:y+2, face)
					this.J_to_sph(2,2,x,y,face) = this.partial_c4(temp, delta)
					this.J_to_cube(2,2,x,y,face) = 2d0*delta/(temp(1) - temp(-1))

				end do
			end do
		end do


		call this.hem_of_face(this.J_to_sph(1, 1, :, :, 1:6))
		call this.hem_of_face(this.J_to_sph(1, 2, :, :, 1:6))
		call this.hem_of_face(this.J_to_sph(2, 1, :, :, 1:6))
		call this.hem_of_face(this.J_to_sph(2, 2, :, :, 1:6))

		call this.hem_of_face(this.J_to_cube(1, 1, :, :, 1:6))
		call this.hem_of_face(this.J_to_cube(1, 2, :, :, 1:6))
		call this.hem_of_face(this.J_to_cube(2, 1, :, :, 1:6))
		call this.hem_of_face(this.J_to_cube(2, 2, :, :, 1:6))

		! this.J_to_sph(:, :, :, :, 3) = this.J_to_sph(:, :, :, :, 2);  this.J_to_sph(:, :, :, :, 4) = this.J_to_sph(:, :, :, :, 2); this.J_to_sph(:, :, :, :, 5) = this.J_to_sph(:, :, :, :, 2)


	end subroutine



	subroutine hem_of_face(this, cubic)
		Class(metric) :: this
		Real(8), intent(inout) :: cubic(1 - this.step:2*this.dim + this.step, 1 - this.step:2*this.dim + this.step, 1:6)
		integer(4) dim, f, l, x, y, face, i, j, k, step

		dim = this.dim;  step = this.step
		f = 1-step; l = 2*dim + step

		do i = 1, step
			do j = 1 - step, 2*dim + step

				cubic(j,2*dim+i,2) = cubic(j,i,6);      cubic(j,1-i,6) = cubic(j,2*dim+1-i,2)
				cubic(2*dim+i,j,2) = cubic(i,j,3);      cubic(1-i,j,3) = cubic(2*dim+1-i,j,2)
				cubic(j,1-i,2) = cubic(j,2*dim+1-i,1);  cubic(j,2*dim+i,1) = cubic(j,i,2)
				cubic(1-i,j,2) = cubic(2*dim+1-i,j,5);  cubic(2*dim+i,j,5) = cubic(i,j,2)

			end do
		end do

		do i = 1, step
			do j = 1 - step, 2*dim + step
			k = 2*dim + 1 -j

				cubic(j,2*dim+i,4) = cubic(k,2*dim+1-i,6);     cubic(k,2*dim+i,6) = cubic(j,2*dim+1-i,4)
				cubic(2*dim+i,j,4) = cubic(i,j,5);             cubic(1-i,j,5) = cubic(2*dim+1-i,j,4)
				cubic(j,1-i,4) = cubic(k,i,1);                 cubic(k,1-i,1) = cubic(j,i,4)
				cubic(1-i,j,4) = cubic(2*dim+1-i,j,3);         cubic(2*dim+i,j,3) = cubic(i,j,4)

			end do
		end do

		do i = 1, step
			do j = 1 - step, 2*dim + step
			k = 2*dim + 1 -j

				cubic(j,2*dim+i,3) = cubic(2*dim+1-i,j,6);     cubic(2*dim+i,j,6) = cubic(j,2*dim+1-i,3)
				cubic(j,1-i,3) = cubic(2*dim+1-i,k,1);         cubic(2*dim+i,k,1) = cubic(j,i,3)

			end do
		end do

		do i = 1, step
			do j = 1 - step, 2*dim + step
			k = 2*dim + 1 -j

				cubic(j,2*dim+i,5) = cubic(i,k,6);     cubic(1-i,k,6) = cubic(j,2*dim+1-i,5)
				cubic(j,1-i,5) = cubic(i,j,1);         cubic(1-i,k,1) = cubic(j,i,5)

			end do
		end do

	end subroutine



	real(8) function partial_c4(this, fun, h)
		Class(metric) :: this
		real(8), intent(in) :: fun(-2:2), h
		real(8) A , B, C, D, E

		A = 2.0/(3.0*h);  B = - 2.0/(3.0*h);  C = - 1.0/(12.0*h);  D = 1.0/(12.0*h);  E = 0.0
		partial_c4 = A*fun(1) + B*fun(-1) + C*fun(2) + D*fun(-2) +  E*fun(0)

	end function



end module