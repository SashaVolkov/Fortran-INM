module geometry

	implicit none

CONTAINS


	subroutine distance_sphere(radius, angl1, angl2, dist)
		real(8), intent(in) :: angl1(1:2),angl2(1:2), radius
		real(8), intent(out) :: dist
		real(8) r1(1:3), r2(1:3)

		dist = radius * dacos(dsin(angl1(1))*dsin(angl2(1)) + dcos(angl1(1))*dcos(angl2(1))*dcos(angl1(2) - angl2(2)))

! 		Кратчайшее расстояние между двумя точками на земной поверхности (если принять ее за сферу) определяется зависимостью:
! cos(d) = sin(φА)·sin(φB) + cos(φА)·cos(φB)·cos(λА − λB),
 ! где φА и φB — широты, λА, λB — долготы данных пунктов, d — расстояние между пунктами, измеряемое в радианах длиной дуги большого круга земного шара. 

! Расстояние между пунктами, измеряемое в километрах, определяется по формуле:
! L = d·R,
! где R = 6371 км — средний радиус земного шара.

	end subroutine


end module geometry
