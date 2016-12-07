module sphere_geometry

implicit none

	Private
	Public :: geometry

	Type geometry

		CONTAINS
		Procedure, Public :: angle => angle_sphere
		Procedure, Public :: triangle => spherical_triangle
	End Type

CONTAINS


	real(8) function angle_sphere(this, latlon1, latlon2)
		Class(geometry) :: this
		real(8), intent(in) :: latlon1(1:2),latlon2(1:2)

		angle_sphere = dacos(dsin(latlon1(1))*dsin(latlon2(1)) + dcos(latlon1(1))*dcos(latlon2(1))*dcos(latlon1(2) - latlon2(2)))

	end function



		subroutine spherical_triangle(this, latlon1, latlon2, latlon3, area, d)

		Class(geometry) :: this
		real(8), intent(in) :: latlon1(1:2), latlon2(1:2), latlon3(1:2)  ! points_latlon
		real(8), intent(out) :: area, d
		real(8) a, b, c  ! angles between radiuses
		real(8) alpha, beta, gamma, eps
		real(8), parameter :: pi = 314159265358979323846d-20


		a = this.angle(latlon2, latlon3)
		b = this.angle(latlon1, latlon3)
		c = this.angle(latlon1, latlon2)

		alpha = dacos( ( dcos(a) - dcos(b)*dcos(c) )/( dsin(b)*dsin(c) ) )
		beta = dacos( ( dcos(b) - dcos(a)*dcos(c) )/( dsin(a)*dsin(c) ) )
		gamma = dacos( ( dcos(c) - dcos(b)*dcos(a) )/( dsin(b)*dsin(a) ) )

		! triangle_angles = (/alpha, beta, gamma/)

		eps = alpha + beta + gamma - pi
		area = eps

		d = c*dsin(alpha)*dsin(beta)/dsin(alpha + beta)

	end subroutine


end module sphere_geometry
