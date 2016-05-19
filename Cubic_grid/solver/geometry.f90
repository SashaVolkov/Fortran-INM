module sphere_geometry

implicit none

	Private
	Public :: geometry

	Type geometry

	real(8) radius, pi

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: dist => distance_sphere
		Procedure, Public :: angle => angle_sphere
		Procedure, Public :: triangle => spherical_triangle
		Procedure, Public :: triangle_area => spherical_triangle_area
	End Type

CONTAINS


	subroutine init(this, radius, pi)
		Class(geometry) :: this
		real(8), intent(in) :: radius, pi

		this.pi = pi;  this.radius = radius

	end subroutine



	subroutine distance_sphere(this, angl1, angl2, dist)
		Class(geometry) :: this
		real(8), intent(in) :: angl1(1:2), angl2(1:2)
		real(8), intent(out) :: dist

		dist = this.radius * dacos(dsin(angl1(1))*dsin(angl2(1)) + dcos(angl1(1))*dcos(angl2(1))*dcos(angl1(2) - angl2(2)))

	end subroutine




	subroutine angle_sphere(this, angl1, angl2, angle)
		Class(geometry) :: this
		real(8), intent(in) :: angl1(1:2),angl2(1:2)
		real(8), intent(out) :: angle

		angle = dacos(dsin(angl1(1))*dsin(angl2(1)) + dcos(angl1(1))*dcos(angl2(1))*dcos(angl1(2) - angl2(2)))

	end subroutine



	!!!!Wiki article "Решение треугольников"
	subroutine spherical_triangle(this, angl1, angl2, angl3, alpha, beta, gamma)

		Class(geometry) :: this
		real(8), intent(in) :: angl1(1:2), angl2(1:2), angl3(1:2)  ! points_latlon
		real(8), intent(out) :: alpha, beta, gamma  ! angles of spherical triangle
		real(8) a, b, c  ! angles between radiuses

		call this.angle(angl1, angl2, a)
		call this.angle(angl2, angl3, b)
		call this.angle(angl1, angl3, c)

		alpha = dacos( ( dcos(a) - dcos(b)*dcos(c) )/( dsin(b)*dsin(c) ) )
		beta = dacos( ( dcos(b) - dcos(a)*dcos(c) )/( dsin(a)*dsin(c) ) )
		gamma = dacos( ( dcos(c) - dcos(b)*dcos(a) )/( dsin(b)*dsin(a) ) )

	end subroutine



		subroutine spherical_triangle_area(this, angl1, angl2, angl3, area)

		Class(geometry) :: this
		real(8), intent(in) :: angl1(1:2), angl2(1:2), angl3(1:2)  ! points_latlon
		real(8), intent(out) :: area
		real(8) a, b, c  ! angles between radiuses
		real(8) alpha, beta, gamma, eps

		call this.angle(angl1, angl2, a)
		call this.angle(angl2, angl3, b)
		call this.angle(angl1, angl3, c)

		alpha = dacos( ( dcos(a) - dcos(b)*dcos(c) )/( dsin(b)*dsin(c) ) )
		beta = dacos( ( dcos(b) - dcos(a)*dcos(c) )/( dsin(a)*dsin(c) ) )
		gamma = dacos( ( dcos(c) - dcos(b)*dcos(a) )/( dsin(b)*dsin(a) ) )

		eps = alpha + beta + gamma - this.pi
		area = this.radius * this.radius * eps

	end subroutine


end module sphere_geometry
