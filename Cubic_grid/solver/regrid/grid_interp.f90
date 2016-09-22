module grid_interp

	use scan_print, Only: printer
	use sphere_geometry, Only: geometry

implicit none

	Private
	Public :: interp

	Type interp

		Real(8), Allocatable :: latlon_c_off(:, :, :, :)
		Real(8), Allocatable :: latlon_c_to(:, :, :)
		Real(8), Allocatable :: weight(:, :, :)
		Real(8), Allocatable :: surface_off(:, :, :)
		Real(8), Allocatable :: surface_to(:, :)
		integer(4), Allocatable :: indexes_xyface(:, :, :, :)
		integer(4), Allocatable :: closest_xyface(:, :, :)

		integer(4) dim, lon_max, lat_max

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Public :: weight_find => weight_find
		Procedure, Public :: hem_of_face => hem_of_face
		Procedure, Public :: nearest_point_search => nearest_point_search
		Procedure, Public :: interpolate => interpolate
	End Type


CONTAINS


	subroutine init(this, dim)

		Class(interp) :: this
		integer(4), intent(in) :: dim

		this.dim = dim;  this.lon_max = 180;  this.lat_max = 90

		call this.alloc()

	end subroutine



	subroutine alloc(this)
		Class(interp) :: this
		integer(4) dim, f, l, lon, lat

		dim = this.dim;  f = 1 - 1; l = 2*dim + 1
		lon = this.lon_max; lat = this.lat_max

		Allocate(this.latlon_c_off(1:2, f:l, f:l, 1:6))
		Allocate(this.surface_off(f:l, f:l, 1:6))
		Allocate(this.surface_to(-lat:lat, -lon:lon))
		Allocate(this.weight(1:4, -lat:lat, -lon:lon))
		Allocate(this.indexes_xyface(1:3, 1:4, -lat:lat, -lon:lon))
		Allocate(this.closest_xyface(1:3, -lat:lat, -lon:lon))

	end subroutine



	subroutine deinit(this)
		Class(interp) :: this
		if (Allocated(this.latlon_c_off)) Deallocate(this.latlon_c_off)
		if (Allocated(this.surface_off)) Deallocate(this.surface_off)
		if (Allocated(this.surface_to)) Deallocate(this.surface_to)
		if (Allocated(this.weight)) Deallocate(this.weight)
		if (Allocated(this.indexes_xyface)) Deallocate(this.indexes_xyface)
	end subroutine



	subroutine weight_find(this, g)
		Class(interp) :: this
		Class(geometry) :: g

		call this.nearest_point_search(g)
		call this.hem_of_face(this.latlon_c_off(1, :, :, :))
		call this.hem_of_face(this.latlon_c_off(2, :, :, :))

	end subroutine


	subroutine nearest_point_search(this, g)
		Class(interp) :: this
		Class(geometry) :: g
		integer(4) dim, f, l, lon, lat, x, y, face(3), i
		Real(8) angle, latlon(1:2), min
		real(8), parameter :: pi = 314159265358979323846d-20

		dim = this.dim
		! lon = this.lon_max; lat = this.lat_max

		do lon = -this.lon_max, this.lon_max
			do lat = -this.lat_max, this.lat_max
			latlon(1) = lat*pi/180d0;  latlon(2) = lon*pi/180d0; min = 10000.0

			if(abs(lon) < 45) then
				face(1:3) = (/1,2,6/)
			else if(abs(lon) > 135) then
				face(1:3) = (/1,4,6/)
			else if(lon < -45)then
				face(1:3) = (/1,5,6/)
			else
				face(1:3) = (/1,3,6/)
			end if

					do x = 1, 2*dim
						do y = 1, 2*dim
							if(abs(lat) < 40)then
angle = g.angle(this.latlon_c_off(1:2, x, y, face(2)), latlon)
if(angle < min) then
	min = angle;  this.closest_xyface(1:3, lat, lon) = (/x,y,face(2)/)
end if

							else if(lat < -50)then
angle = g.angle(this.latlon_c_off(1:2, x, y, face(1)), latlon)
if(angle < min) then
	min = angle;  this.closest_xyface(1:3, lat, lon) = (/x,y,face(1)/)
end if

							else if(lat > 50)then
angle = g.angle(this.latlon_c_off(1:2, x, y, face(3)), latlon)
if(angle < min) then
	min = angle;  this.closest_xyface(1:3, lat, lon) = (/x,y,face(3)/)
end if

							else
do i = 1, 3
	angle = g.angle(this.latlon_c_off(1:2, x, y, face(i)), latlon)
	if(angle < min) then
		min = angle;  this.closest_xyface(1:3, lat, lon) = (/x,y,face(i)/)
	end if
end do

							end if
					end do
				end do
			end do
		end do
	end subroutine


	subroutine interpolate(this)
		Class(interp) :: this
		integer(4) dim, f, l, lon, lat, x, y, face, i

		! surf_to = sum(w*surf_off)
		call this.hem_of_face(this.surface_off)

		do lon = -this.lon_max, this.lon_max
			do lat = -this.lat_max, this.lat_max
				this.surface_to(lat, lon) = 0d0
				do i = 1, 4
					x = this.indexes_xyface(1, i, lat, lon)
					y = this.indexes_xyface(2, i, lat, lon)
					face = this.indexes_xyface(3, i, lat, lon)
					this.surface_to(lat, lon) = this.surface_to(lat, lon) + this.surface_off(x,y,face)*this.weight(i, lat, lon)
				end do
			end do
		end do

	end subroutine



	subroutine hem_of_face(this, cubic)
		Class(interp) :: this
		Real(8), intent(inout) :: cubic(0:2*this.dim+1, 0:2*this.dim+1, 1:6)
		integer(4) dim, f, l, x, y, face, i

		dim = this.dim;

		cubic(2*dim+1,:,2) = cubic(1,:,3);  cubic(0,:,3) = cubic(2*dim,:,2)
		cubic(2*dim+1,:,3) = cubic(1,:,4);  cubic(0,:,4) = cubic(2*dim,:,3)
		cubic(2*dim+1,:,4) = cubic(1,:,5);  cubic(0,:,5) = cubic(2*dim,:,4)
		cubic(2*dim+1,:,5) = cubic(1,:,2);  cubic(0,:,2) = cubic(2*dim,:,5)

		cubic(:,2*dim+1,2) = cubic(:,1,6);  cubic(:,0,6) = cubic(:,2*dim,2)
		cubic(:,2*dim+1,3) = cubic(2*dim,:,6);  cubic(2*dim+1,:,6) = cubic(:,2*dim,3)
		cubic(:,2*dim+1,4) = cubic(2*dim:1,2*dim,6);  cubic(2*dim:1,2*dim+1,6) = cubic(:,2*dim,4)
		do i = 1, 2*dim
			cubic(i,2*dim+1,5) = cubic(1,2*dim + 1 - i,6);  cubic(0,2*dim + 1 - i,6) = cubic(i,2*dim,5)
		end do

		cubic(:,0,2) = cubic(:,2*dim,1);  cubic(:,2*dim+1,1) = cubic(:,1,2)
		cubic(:,0,3) = cubic(2*dim,2*dim:1,1);  cubic(2*dim+1,2*dim:1,1) = cubic(:,1,3)
		cubic(:,0,4) = cubic(2*dim:1,1,1);  cubic(2*dim:1,0,1) = cubic(:,1,4)
		cubic(:,0,5) = cubic(1,:,1);  cubic(0,:,1) = cubic(:,1,5)

	end subroutine



end module