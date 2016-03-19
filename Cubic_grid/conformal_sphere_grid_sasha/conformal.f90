!!!!To know more look through Rancic Parser Mesinger article 1996 Apendix A & B!!!!!!! 

module conformal

IMPLICIT NONE

	Private
	Public :: conf

	Type conf
		CONTAINS
		Procedure :: conformal_z_w => conformal_z_w
	End Type

CONTAINS

subroutine conformal_z_w(this, z, w)

	Class(conf) :: this

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

	BIG_Z = z ** 4

	BIG_W = 0
	do i=1,30
		 BIG_W = BIG_W + A_c(i)*(BIG_Z**i)
	end do
	w = BIG_W ** (3333333333333333d-16)

end subroutine conformal_z_w


subroutine conformal_w_z(this, w,z)

	Class(conf) :: this
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

  BIG_W = w ** 3
  BIG_Z = 0
  do i=1,30
     BIG_Z = BIG_Z + B_c(i)*(BIG_W**i)
  end do
  z = BIG_Z ** (25d-2)

end subroutine conformal_w_z


end module conformal