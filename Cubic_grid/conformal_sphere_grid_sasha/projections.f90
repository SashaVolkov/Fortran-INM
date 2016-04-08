module projections

implicit none
Private
Public :: projection

Type projection
	CONTAINS
	Procedure :: inverse => inverse_stereo_projection
	Procedure :: conformal_z_w => conformal_z_w
	Procedure :: conformal_w_z => conformal_w_z
	Procedure :: stereographic_cube_to_sphere => cube2sphere
	Procedure :: direct => direct_stereo_projection
End Type


CONTAINS

	subroutine inverse_stereo_projection(this, x_plane, y_plane,r_sphere, x, y, z, status)
	! reverse central projection from south pole of sphere
	! with centre in (0,0,0) on tangent plane in north pole (0,0,R)
	! return: 
	!  0 if ok
	!  1 if error, r_sphere < 0
		real(8), intent(in) :: x_plane, y_plane, r_sphere
		real(8), intent(out) :: x, y, z
		integer(4), intent(out) :: status
		Class(projection) :: this

		real(8) r_coeff

		r_coeff = r_sphere / (2d0 + x_plane**2 +y_plane**2)

		if (r_sphere < 0d0) then
			write(*,*) "Warning, radius of sphere is less than 0:"
			write(*,*) "radius =", r_sphere
			status = 1
		else
			x = 2*sqrt(2d0) * x_plane * r_coeff
			y = 2*sqrt(2d0) * y_plane * r_coeff
			z = (2-x_plane**2-y_plane**2) * r_coeff
			status = 0
		end if
	end subroutine



	subroutine conformal_z_w(this, z, w)
	! conformal projection z -> w 

	!!!!To know more look through Rancic Parser Mesinger article 1996 Apendix A & B!!!!!!! 
	!!!!Q.J.R. Meteorol. Soc. 1996, 122, pp. 959-982!!!!!

		Class(projection) :: this

		real(8), dimension(:), parameter :: A_c(1:30) = &
				 (/       147713062600964d-14, &
			-38183510510174d-14, &
			-05573058001191d-14, &
			-00895883606818d-14, &
			-00791315785221d-14, & ! 5
			-00486625437708d-14, &
			-00329251751279d-14, &
			-00235481488325d-14, &
			-00175870527475d-14, &
			-00135681133278d-14, & !10
			-00107459847699d-14, &
			-00086944475948d-14, &
			-00071607115121d-14, &
			-00059867100093d-14, &
			-00050699063238d-14, & !15
			-00043415191279d-14, &
			-00037541003286d-14, &
			-00032741060100d-14, &
			-00028773091482d-14, &
			-00025458777519d-14, & !20
			-00022664642371d-14, & 
			-00020289261022d-14, &
			-00018254510830d-14, &
			-00016499474461d-14, &
			-00014976117168d-14, & !25
			-00013646173946d-14, &
			-00012478875823d-14, &
			-00011449267279d-14, &
			-00010536946150d-14, &
			-00009725109376d-14 /)  !!Page 981 Table B1

		complex*16 w, z, BIG_W, BIG_Z
		integer i

		BIG_Z = z ** 4  !!Page 977 eq. A.2

		BIG_W = 0
		do i=1,30
			 BIG_W = BIG_W + A_c(i)*(BIG_Z**i)  !!Page 977 eq. A.3a
		end do
		w = BIG_W ** (3333333333333333d-16)

	end subroutine conformal_z_w


	subroutine conformal_w_z(this, w,z)
	! conformal projection w -> z

		Class(projection) :: this
	  complex*8 w,z, BIG_W,BIG_Z
	  integer i
	  real(8), dimension(:), parameter:: &
	       B_c(30)=(/.67698819751739, &
	    .11847293456564, &
	    .05317178134668, &
	    .02965810434052, &
	    .01912447304028, &
	    .01342565621117, &
	    .00998873323180, &
	    .00774868996406, &
	    .00620346979888, &
	    .00509010874883, &
	    .00425981184328, &
	    .00362308956077, &
	    .00312341468940, &
	    .00272360948942, &
	    .00239838086555, &
	    .00213001905118, &
	    .00190581316131, &
	    .00171644156404, &
	    .00155493768255, &
	    .00141600715207, &
	    .00129556597765, &
	    .00119042140226, &
	    .00109804711790, &
	    .00101642216628, &
	    .00094391366522, &
	    .00087919021224, &
	    .00082115710311, &
	    .00076890728775, &
	    .00072168382969, &
	    .00067885087750  /)  !!Page 981 Table B1

	  BIG_W = w ** 3  !!Page 978 (iii)
	  BIG_Z = 0
	  do i=1,30
	     BIG_Z = BIG_Z + B_c(i)*(BIG_W**i)  !!Page 978 (iv)
	  end do
	  z = BIG_Z ** (25d-2)  !!Page 978 (v)

	end subroutine conformal_w_z


	subroutine cube2sphere(this, x, y, z, x_face, y_face, r_sphere, face_index, status)
	 ! stereographic projection between points on face of cubed sphere on real sphere
	 ! xi_face - coordinate on number's face, in [-1, 1]
	 ! r_sphere - radius of inscribed sphere
	 ! x, y, z - cartesian coordinates of point

		real(8), intent(in) :: x_face, y_face, r_sphere
		integer(4), intent(in) :: face_index

		real(8), intent(out) :: x, y, z
		integer(4), intent(out) ::  status

		real(8) BIG_R
		Class(projection) :: this

		status = 0
		BIG_R = r_sphere / sqrt(1 + x_face**2 + y_face**2) ! Baiburin p.5 (3.1)

		select case (face_index)
			 case (1)
					x = y_face * BIG_R
					y = x_face * BIG_R		! Baiburin (3.2)
					z = - BIG_R
			 case (2)
					x = BIG_R
					y = x_face * BIG_R		! Baiburin (3.3)
					z = y_face * BIG_R
			 case (3)
					x = - x_face * BIG_R
					y = BIG_R							! Baiburin (3.4)
					z = y_face * BIG_R
			 case (4)
					x = - BIG_R
					y = - x_face * BIG_R	! Baiburin (3.5)
					z = y_face * BIG_R
			 case (5)
					x = x_face * BIG_R
					y = - BIG_R						! Baiburin (3.6)
					z = y_face * BIG_R
			 case (6)
					x = -y_face * BIG_R
					y = x_face * BIG_R		! Baiburin (3.7)
					z = BIG_R 
			 case default
					status = 1
		end select

	end subroutine


	subroutine direct_stereo_projection(this, x_plane, y_plane,r_sphere, x, y, z)
	! central projection from south pole of sphere
	! with centre in (0,0,0) on tangent plane in north pole (0,0,R)
		real(8), intent(out) :: x_plane, y_plane, r_sphere
		real(8), intent(in) :: x, y, z
		Class(projection) :: this

		r_sphere = sqrt(x**2+y**2+z**2)
		x_plane = sqrt(2d0)*x/(r_sphere+z)
		y_plane = sqrt(2d0)*y/(r_sphere+z)
	end subroutine

end module projections