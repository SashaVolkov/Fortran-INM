module grid_interp

	use scan_print, Only: printer
	use sphere_geometry, Only: geometry
	use omp_lib

implicit none

	Private
	Public :: interp

	Type interp

		Real(8), Allocatable :: latlon_c_off(:, :, :, :)
		Real(8), Allocatable :: latlon_c_to(:, :, :)
		Real(8), Allocatable :: weight(:,:, :, :)
		Real(8), Allocatable :: surface_off(:, :, :, :)
		Real(8), Allocatable :: surface_to(:, :, :)
		integer(4), Allocatable :: indexes_xyface(:, :, :, :)
		integer(4), Allocatable :: closest_xyface(:, :, :)

		integer(4) dim, lon_max, lat_max, step

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Public :: weight_find => weight_find
		Procedure, Public :: hem_of_face => hem_of_face
		Procedure, Public :: nearest_point_search => nearest_point_search
		Procedure, Public :: cell_search => cell_search
		Procedure, Public :: interpolate => interpolate
	End Type


CONTAINS


	subroutine init(this, dim)

		Class(interp) :: this
		integer(4), intent(in) :: dim

		this.dim = dim;  this.lon_max = 180;  this.lat_max = 90
		this.step = 2

		call this.alloc()

	end subroutine



	subroutine alloc(this)
		Class(interp) :: this
		integer(4) dim, f, l, lon, lat

		dim = this.dim;  f = 1 - this.step; l = 2*dim + this.step
		lon = this.lon_max; lat = this.lat_max

		Allocate(this.latlon_c_off(1:2, f:l, f:l, 1:6))
		Allocate(this.surface_off(f:l, f:l, 1:6, 2))
		Allocate(this.surface_to(-lon:lon, -lat:lat, 2))
		Allocate(this.weight(4, 4, -lat:lat, -lon:lon))
		Allocate(this.indexes_xyface(1:3, 1:4, -lat:lat, -lon:lon))
		Allocate(this.closest_xyface(1:3, -lon:lon, -lat:lat))

	end subroutine



	subroutine deinit(this)
		Class(interp) :: this
		if (Allocated(this.latlon_c_off)) Deallocate(this.latlon_c_off)
		if (Allocated(this.surface_off)) Deallocate(this.surface_off)
		if (Allocated(this.surface_to)) Deallocate(this.surface_to)
		if (Allocated(this.weight)) Deallocate(this.weight)
		if (Allocated(this.indexes_xyface)) Deallocate(this.indexes_xyface)
	end subroutine



	subroutine weight_find(this, g, printer_scaner)
		Class(interp) :: this
		Class(printer) :: printer_scaner
		Class(geometry) :: g
		integer(4) lon, lat, x(4), y(4), face(4), i, j, dim
		Real(8) :: S(4), Big_S, latlon(2), latlon1(2), latlon2(2), d(4), ang(4), d_xy(2,4), lagr(2,4)
		real(8), parameter :: pi = 314159265358979323846d-20

		call this.hem_of_face(this.latlon_c_off(1, :, :, :))
		call this.hem_of_face(this.latlon_c_off(2, :, :, :))
		if(printer_scaner.point_find == 0) then
			call printer_scaner.scan_point(this.closest_xyface)
		else if(printer_scaner.point_find == 1) then
			call this.nearest_point_search(g)
			call printer_scaner.print_point(this.closest_xyface)
		end if

		call this.cell_search(g)

		dim = this.dim

		do lon = -this.lon_max, this.lon_max
			do lat = -this.lat_max, this.lat_max

				latlon(1) = lat*pi/180d0;  latlon(2) = lon*pi/180d0

				face(:) = this.indexes_xyface(3, :, lat, lon)

				x(:) = this.indexes_xyface(1, :, lat, lon)
				y(:) = this.indexes_xyface(2, :, lat, lon)
				

				latlon1(:) = this.latlon_c_off(:, x(1), y(1), face(1))
				latlon2(:) = this.latlon_c_off(:, x(2), y(2), face(2))
				call g.triangle(latlon1, latlon2, latlon, S(1), d(1))
				ang(1) = g.angle(latlon1, latlon2)

				latlon1(:) = this.latlon_c_off(:, x(2), y(2), face(2))
				latlon2(:) = this.latlon_c_off(:, x(3), y(3), face(3))
				call g.triangle(latlon1, latlon2, latlon, S(2), d(2))
				ang(2) = g.angle(latlon1, latlon2)

				latlon1(:) = this.latlon_c_off(:, x(3), y(3), face(3))
				latlon2(:) = this.latlon_c_off(:, x(4), y(4), face(4))
				call g.triangle(latlon1, latlon2, latlon, S(3), d(3))
				ang(3) = g.angle(latlon1, latlon2)

				latlon1(:) = this.latlon_c_off(:, x(4), y(4), face(4))
				latlon2(:) = this.latlon_c_off(:, x(1), y(1), face(1))
				call g.triangle(latlon1, latlon2, latlon, S(4), d(4))
				ang(4) = g.angle(latlon1, latlon2)

				! Big_S = (S(1) + S(3))*(S(2) + S(4))

				! this.weight(1, lat, lon) = S(2)*S(3)/Big_S
				! this.weight(2, lat, lon) = S(4)*S(3)/Big_S
				! this.weight(3, lat, lon) = S(1)*S(4)/Big_S
				! this.weight(4, lat, lon) = S(2)*S(1)/Big_S

				d_xy(1,1) = ang(1) + d(4);  d_xy(1,2) = d(4);  d_xy(1,3) = d(2);  d_xy(1,4) = ang(1) + d(2)
				d_xy(2,1) = ang(4) + d(3);  d_xy(2,2) = d(3);  d_xy(2,3) = d(1);  d_xy(2,4) = ang(4) + d(1)

				Big_S = (d(1) + d(3))*(d(2) + d(4))

				this.weight(2,3, lat, lon) = d(2)*d(3)/Big_S
				this.weight(3,3, lat, lon) = d(4)*d(3)/Big_S
				this.weight(3,2, lat, lon) = d(1)*d(4)/Big_S
				this.weight(2,2, lat, lon) = d(2)*d(1)/Big_S

				if (x(4) > 1 .and. y(4) > 1 .and. y(2) < 2*dim .and. x(2) < 2*dim) then
					lagr(:,1) =-(d_xy(:,2)*d_xy(:,3)*d_xy(:,4))/((d_xy(:,1) - d_xy(:,2))*(d_xy(:,1) + d_xy(:,3))*(d_xy(:,1) + d_xy(:,4)))
					lagr(:,2) = (d_xy(:,1)*d_xy(:,3)*d_xy(:,4))/((d_xy(:,1) - d_xy(:,2))*(d_xy(:,2) + d_xy(:,3))*(d_xy(:,2) + d_xy(:,4)))
					lagr(:,3) = (d_xy(:,1)*d_xy(:,2)*d_xy(:,4))/((d_xy(:,1) + d_xy(:,3))*(d_xy(:,2) + d_xy(:,3))*(d_xy(:,4) - d_xy(:,3)))
					lagr(:,4) =-(d_xy(:,1)*d_xy(:,2)*d_xy(:,3))/((d_xy(:,1) + d_xy(:,4))*(d_xy(:,2) + d_xy(:,4))*(d_xy(:,4) - d_xy(:,3)))

					do i = 1, 4
						do j = 1, 4
							this.weight(i,j, lat, lon) = lagr(1, i)*lagr(2, j)
						end do
					end do
				end if
			end do
		end do
		! print *, this.weight(:,:,20,20)
		! print *, sum(this.weight(:,:,20,20))
		! print *, sum(lagr(1, :))

	end subroutine


	subroutine nearest_point_search(this, g)
		Class(interp) :: this
		Class(geometry) :: g
		integer(4) dim, f, l, lon, lat, x, y, face(3), i
		Real(8) angle, latlon(1:2), min
		real(8), parameter :: pi = 314159265358979323846d-20

		dim = this.dim
		! lon = this.lon_max; lat = this.lat_max

		!$OMP PARALLEL PRIVATE(i, y, x, lon, lat, angle, latlon, min, face)
		!$OMP DO

		do lon = -this.lon_max, this.lon_max
			do lat = -this.lat_max, this.lat_max
			latlon(1) = lat*pi/180d0;  latlon(2) = lon*pi/180d0; min = 10000.0

					do x = 1, 2*dim
						do y = 1, 2*dim
							if(lat < -50)then
								angle = g.angle(this.latlon_c_off(1:2, x, y, 1), latlon)
								if(angle < min) then
									min = angle;  this.closest_xyface(1:3, lon, lat) = (/x,y,1/)
								end if
							else if(lat > 50)then
								angle = g.angle(this.latlon_c_off(1:2, x, y, 6), latlon)
								if(angle < min) then
									min = angle;  this.closest_xyface(1:3, lon, lat) = (/x,y,6/)
								end if
							else

					if(abs(lon) <= 45) then
						face(1:3) = (/1,2,6/)
					else if(abs(lon) >= 135) then
						face(1:3) = (/1,4,6/)
					else if(lon < -45 .and. lon > -135)then
						face(1:3) = (/1,5,6/)
					else if(lon > 45 .and. lon < 135)then
						face(1:3) = (/1,3,6/)
					end if

								do i = 1, 3
									if( abs(this.latlon_c_off(2, x, y, face(i)) - latlon(2)) < 0.5 .and. abs(this.latlon_c_off(1, x, y, face(i)) - latlon(1)) < 0.5 ) then
										angle = g.angle(this.latlon_c_off(1:2, x, y, face(i)), latlon)
										if(angle < min .and. angle >= 0d0) then
											min = angle;  this.closest_xyface(1:3, lon, lat) = (/x,y,face(i)/)
										end if
									end if
								end do
							end if
					end do
				end do
			end do
		end do

		!$OMP END DO
		!$OMP END PARALLEL

	end subroutine


	subroutine interpolate(this)
		Class(interp) :: this
		integer(4) dim, f, l, lon, lat, x(4), y(4), face, i, j
		Real(8) :: weight(4, 4)

		! surf_to = sum(w*surf_off)
		call this.hem_of_face(this.surface_off(:,:,:,1))
		call this.hem_of_face(this.surface_off(:,:,:,2))

		!$OMP PARALLEL PRIVATE(i, lon, lat, face)
		!$OMP DO

		do lon = -this.lon_max, this.lon_max
			do lat = -this.lat_max, this.lat_max
				this.surface_to(lon, lat, :) = 0d0
					x(:) = this.indexes_xyface(1, :, lat, lon)
					y(:) = this.indexes_xyface(2, :, lat, lon)
					face = this.closest_xyface(3, lon, lat)
					weight(:,:) = this.weight(:,:, lat, lon)

this.surface_to(lon, lat, :) = this.surface_off(x(1),y(1),face, :)*weight(2,3) + this.surface_off(x(2),y(2),face, :)*weight(3,3) + this.surface_off(x(3),y(3),face, :)*weight(3,2) + this.surface_off(x(4),y(4),face, :)*weight(2,2) + &
this.surface_off(x(1)-1,y(1)+1,face, :)*weight(1,4) + this.surface_off(x(1),y(1)+1,face, :)*weight(2,4) + this.surface_off(x(2),y(1)+1,face, :)*weight(3,4) + this.surface_off(x(2)+1,y(1)+1,face, :)*weight(4,4) + &
this.surface_off(x(1)-1,y(1),face, :)*weight(1,3) + this.surface_off(x(1)-1,y(4),face, :)*weight(1,2) + this.surface_off(x(2)+1,y(2),face, :)*weight(4,3) + this.surface_off(x(3)+1,y(3),face, :)*weight(4,2) + &
this.surface_off(x(4)-1,y(4)-1,face, :)*weight(1,1) + this.surface_off(x(4),y(4)-1,face, :)*weight(2,1) + this.surface_off(x(3),y(3)-1,face, :)*weight(3,1) + this.surface_off(x(3)+1,y(3)-1,face, :)*weight(4,1)

			end do
		end do

		!$OMP END DO
		!$OMP END PARALLEL

	end subroutine



	subroutine hem_of_face(this, cubic)
		Class(interp) :: this
		Real(8), intent(inout) :: cubic(1-this.step:2*this.dim+this.step, 1-this.step:2*this.dim+this.step, 1:6)
		integer(4) dim, f, l, x, y, face, i, j, k, step

		dim = this.dim;  step = this.step

		do i = 1, step
			do j = 1, 2*dim

				cubic(j,2*dim+i,2) = cubic(j,i,6);      cubic(j,1-i,6) = cubic(j,2*dim+1-i,2)
				cubic(2*dim+i,j,2) = cubic(i,j,3);      cubic(1-i,j,3) = cubic(2*dim+1-i,j,2)
				cubic(j,1-i,2) = cubic(j,2*dim+1-i,1);  cubic(j,2*dim+i,1) = cubic(j,i,2)
				cubic(1-i,j,2) = cubic(2*dim+1-i,j,5);  cubic(2*dim+i,j,5) = cubic(i,j,2)

			end do
		end do

		do i = 1, step
			do j = 1, 2*dim
			k = 2*dim + 1 -j

				cubic(j,2*dim+i,4) = cubic(k,2*dim+1-i,6);     cubic(k,2*dim+i,6) = cubic(j,2*dim+1-i,4)
				cubic(2*dim+i,j,4) = cubic(i,j,5);             cubic(1-i,j,5) = cubic(2*dim+1-i,j,4)
				cubic(j,1-i,4) = cubic(k,i,1);                 cubic(k,1-i,1) = cubic(j,i,4)
				cubic(1-i,j,4) = cubic(2*dim+1-i,j,3);         cubic(2*dim+i,j,3) = cubic(i,j,4)

			end do
		end do

		do i = 1, step
			do j = 1, 2*dim
			k = 2*dim + 1 -j

				cubic(j,2*dim+i,3) = cubic(2*dim+1-i,j,6);     cubic(2*dim+i,j,6) = cubic(j,2*dim+1-i,3)
				cubic(j,1-i,3) = cubic(2*dim+1-i,k,1);         cubic(2*dim+i,k,1) = cubic(j,i,3)

			end do
		end do

		do i = 1, step
			do j = 1, 2*dim
			k = 2*dim + 1 -j

				cubic(j,2*dim+i,5) = cubic(i,k,6);     cubic(1-i,k,6) = cubic(j,2*dim+1-i,5)
				cubic(j,1-i,5) = cubic(i,j,1);         cubic(1-i,j,1) = cubic(j,i,5)

			end do
		end do

		do i = 1, step
			do j = 1, step
			cubic(1-i, 1-j, :) = cubic(2*step+1-j, i, :)
			cubic(1-i, 2*dim + j, :) = cubic(1-j, 2*dim+1-i, :)

			cubic(2*dim + i, 1-j, :) = cubic(2*dim+j, i, :)
			cubic(2*dim + i, 2*dim + j, :) = cubic(2*dim+j, 2*dim+1-i, :)
			end do
		end do

	end subroutine



	subroutine cell_search(this, g)
		Class(interp) :: this
		Class(geometry) :: g
		integer(4) dim, f, l, lat, lon, x, y, face, i, x_cell(4), y_cell(4), point
		Real(8) angle, latlon(1:2), min
		real(8), parameter :: pi = 314159265358979323846d-20
		dim = this.dim

		do lon = -this.lon_max, this.lon_max
			do lat = -this.lat_max, this.lat_max

				latlon(1) = lat*pi/180d0;  latlon(2) = lon*pi/180d0; min = 10000.0

				x = this.closest_xyface(1, lon, lat)
				y = this.closest_xyface(2, lon, lat)
				face = this.closest_xyface(3, lon, lat)

				if(x == 1 .and. y == 1) y = 2
				if(x == 1 .and. y == 2*dim) x = 2
				if(x == 2*dim .and. y == 1) x = 2*dim-1
				if(x == 2*dim .and. y == 2*dim) y = 2*dim-1


				x_cell(1) = x - 1; y_cell(1) = y + 1
				x_cell(2) = x + 1; y_cell(2) = y + 1
				x_cell(3) = x + 1; y_cell(3) = y - 1
				x_cell(4) = x - 1; y_cell(4) = y - 1

				if(x_cell(1) == 0 .and. y_cell(1) == 2*dim+1) then
					x_cell(1) = 1
				else if(x_cell(2) == 2*dim+1 .and. y_cell(2) == 2*dim+1) then
					y_cell(2) = 2*dim
				else if(x_cell(3) == 2*dim+1 .and. y_cell(3) == 0) then
					x_cell(3) = 2*dim
				else if(x_cell(4) == 0 .and. y_cell(4) == 0) then
					y_cell(2) = 1
				end if


					do i = 1, 4
						angle = g.angle(this.latlon_c_off(:, x_cell(i), y_cell(i), face), latlon)
						if(angle < min) then
							min = angle;  point = i
						end if
					end do



				if( point == 1) then
					this.indexes_xyface(1:3, 1, lat, lon) = (/x-1,y+1,face/)
					this.indexes_xyface(1:3, 2, lat, lon) = (/x,y+1,face/)
					this.indexes_xyface(1:3, 3, lat, lon) = (/x,y,face/)
					this.indexes_xyface(1:3, 4, lat, lon) = (/x-1,y,face/)
				else if( point == 2) then
					this.indexes_xyface(1:3, 1, lat, lon) = (/x,y+1,face/)
					this.indexes_xyface(1:3, 2, lat, lon) = (/x+1,y+1,face/)
					this.indexes_xyface(1:3, 3, lat, lon) = (/x+1,y,face/)
					this.indexes_xyface(1:3, 4, lat, lon) = (/x,y,face/)
				else if( point == 3) then
					this.indexes_xyface(1:3, 1, lat, lon) = (/x,y,face/)
					this.indexes_xyface(1:3, 2, lat, lon) = (/x+1,y,face/)
					this.indexes_xyface(1:3, 3, lat, lon) = (/x+1,y-1,face/)
					this.indexes_xyface(1:3, 4, lat, lon) = (/x,y-1,face/)
				else if( point == 4) then
					this.indexes_xyface(1:3, 1, lat, lon) = (/x-1,y,face/)
					this.indexes_xyface(1:3, 2, lat, lon) = (/x,y,face/)
					this.indexes_xyface(1:3, 3, lat, lon) = (/x,y-1,face/)
					this.indexes_xyface(1:3, 4, lat, lon) = (/x-1,y-1,face/)
				end if


			end do
		end do

		! this.latlon_c_off(1:2, x_cell, x_cell, face)
	end subroutine



end module