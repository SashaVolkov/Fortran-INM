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
		! Procedure, Public :: triangle_area => spherical_triangle_area
	End Type

CONTAINS


	subroutine init(this, radius, pi)
		Class(geometry) :: this
		real(8), intent(in) :: radius, pi

		this.pi = pi;  this.radius = radius

	end subroutine



	subroutine distance_sphere(this, latlon1, latlon2, dist)
		Class(geometry) :: this
		real(8), intent(in) :: latlon1(1:2), latlon2(1:2)
		real(8), intent(out) :: dist

		dist = this.radius * dacos(dsin(latlon1(2))*dsin(latlon2(2)) + dcos(latlon1(2))*dcos(latlon2(2))*dcos(latlon1(1) - latlon2(1)))

	end subroutine




	subroutine angle_sphere(this, latlon1, latlon2, angle)
		Class(geometry) :: this
		real(8), intent(in) :: latlon1(1:2),latlon2(1:2)
		real(8), intent(out) :: angle

		angle = dacos(dsin(latlon1(2))*dsin(latlon2(2)) + dcos(latlon1(2))*dcos(latlon2(2))*dcos(latlon1(1) - latlon2(1)))

	end subroutine



		subroutine spherical_triangle(this, latlon1, latlon2, latlon3, area, triangle_angles)

		Class(geometry) :: this
		real(8), intent(in) :: latlon1(1:2), latlon2(1:2), latlon3(1:2)  ! points_latlon
		real(8), intent(out) :: area, triangle_angles(1:3)
		real(8) a, b, c  ! angles between radiuses
		real(8) alpha, beta, gamma, eps

		call this.angle(latlon1, latlon2, a)
		call this.angle(latlon2, latlon3, b)
		call this.angle(latlon1, latlon3, c)

		alpha = dacos( ( dcos(a) - dcos(b)*dcos(c) )/( dsin(b)*dsin(c) ) )
		beta = dacos( ( dcos(b) - dcos(a)*dcos(c) )/( dsin(a)*dsin(c) ) )
		gamma = dacos( ( dcos(c) - dcos(b)*dcos(a) )/( dsin(b)*dsin(a) ) )

		triangle_angles = (/alpha, beta, gamma/)

		eps = alpha + beta + gamma - this.pi
		area = this.radius * this.radius * eps

	end subroutine


end module sphere_geometry
