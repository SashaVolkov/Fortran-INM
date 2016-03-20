!!!!To know more look through Rancic Parser Mesinger article 1996 Apendix A & B!!!!!!! 
!!!!Q.J.R. Meteorol. Soc. 1996, 122, pp. 959-982!!!!!

module mapping

IMPLICIT NONE

	Private
	Public :: mapp

	Type mapp
		CONTAINS
		Procedure :: conformal_z_w => conformal_z_w
		Procedure :: conformal_w_z => conformal_w_z
		Procedure :: cube2sphere => cube2sphere
	End Type



CONTAINS

	subroutine conformal_z_w(this, z, w)
	! mapping z -> w 

		Class(mapp) :: this

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
			-00009725109376d-14 /)
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

		Class(mapp) :: this
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
	    .00067885087750  /)

	  BIG_W = w ** 3  !!Page 977 eq. A.2
	  BIG_Z = 0
	  do i=1,30
	     BIG_Z = BIG_Z + B_c(i)*(BIG_W**i)  !!Page 977 eq. A.3b
	  end do
	  z = BIG_Z ** (25d-2)

	end subroutine conformal_w_z


	subroutine cube2sphere(this, x, y, z, x_edge, y_edge, r_sphere, number_edge, status)
	 ! mapping between points on edge of cubed sphere on real sphere
	 ! xi_edge - coordinate on number's edge, in [-1, 1]
	 ! r_sphere - radius of inscribed sphere
	 ! x, y, z - cartesian coordinates of point

		real(8) x, y, z
		real(8) x_edge, y_edge, r_sphere
		integer number_edge, status

		real(8) rr
		Class(mapp) :: this

		status = 0
		rr = r_sphere / sqrt(1 + x_edge**2 + y_edge**2)

		select case (number_edge)
			 case (1)
					x = y_edge * rr
					y = x_edge * rr
					z = - rr
			 case (2)
					x = rr
					y = x_edge * rr
					z = y_edge * rr
			 case (3)
					x = - x_edge * rr
					y = rr
					z = y_edge * rr
			 case (4)
					x = - rr
					y = - x_edge * rr
					z = y_edge * rr
			 case (5)
					x = x_edge * rr
					y = - rr
					z = y_edge * rr
			 case (6)
					x = -y_edge * rr
					y = x_edge * rr
					z = rr 
			 case default
					status = 1
		end select

	end subroutine


end module mapping