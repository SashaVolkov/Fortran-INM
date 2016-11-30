module diagn


implicit none

	Private
	Public :: diagnostic

	Type diagnostic

		Real(4), Allocatable :: surface_precise(:, :)
		integer(4) lon_max, lat_max
		real(8) convert_time

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Private :: alloc => alloc
		Procedure, Public :: deinit => deinit
		Procedure, Public :: L_norm => L_norm
	End Type


CONTAINS


	subroutine init(this, dim, convert_time, grid_type, rescale)

		Class(diagnostic) :: this
		integer(4), intent(in) :: dim, grid_type, rescale
		real(8), intent(in) :: convert_time

		character(32) istring

		this.lon_max = 180;  this.lat_max = 90
		this.convert_time = convert_time/(644d0*4000d0)

		write(istring, *) 2*dim

		if (grid_type == 1) then
			istring = trim(adjustl(istring))//'/equiang/'
		else if (grid_type == 0) then
			if (rescale == 1) then
				istring = trim(adjustl(istring))//'/tan/'
			else if (rescale == 0) then
				istring = trim(adjustl(istring))//'/simple/'
			else if (rescale == 2) then
				istring = trim(adjustl(istring))//'/exp/'
			end if
		end if

		call this.alloc()

		open(11,file='../datFiles/'//trim(istring)//'L1.dat')
		open(12,file='../datFiles/'//trim(istring)//'L2.dat')
		open(13,file='../datFiles/'//trim(istring)//'C.dat')
		open(14,file='../datFiles/'//trim(istring)//'L_inf_cube.dat', FORM="FORMATTED",STATUS="OLD",ACTION="READ")

	end subroutine



	subroutine alloc(this)
		Class(diagnostic) :: this
		integer(4) lon, lat

		lon = this.lon_max; lat = this.lat_max

		Allocate(this.surface_precise(-lon:lon, -lat:lat))

	end subroutine



	subroutine deinit(this)
		Class(diagnostic) :: this
		if (Allocated(this.surface_precise)) Deallocate(this.surface_precise)
		close(11)
		close(12)
		close(13)
		close(14)
	end subroutine



	subroutine L_norm(this,time, surface_to)
		Class(diagnostic) :: this
		integer(4), intent(in) :: time
		real(8), intent(in) :: surface_to(-this.lon_max:this.lon_max, -this.lat_max:this.lat_max)

		integer(4) lat, lon, id, ier
		real(4) L1, L2, L1_prec, L2_prec, L_inf, L_inf_prec, F1, F1_prec, square, pi

		pi = 314159265358979323846d-20

		do lat = -this.lat_max+1, this.lat_max-1
			square = cos(lat*pi/180.0)
			do lon = -this.lon_max+1, this.lon_max-1

				if (isnan(surface_to(lon,lat))) then
				else
				F1 = this.surface_precise(lon, lat) - surface_to(lon,lat)
				F1_prec = this.surface_precise(lon, lat)

				L1 = abs(F1)*square + L1
				L2 = F1*F1*square + L2

				L1_prec = abs(F1_prec)*square + L1_prec
				L2_prec = F1_prec*F1_prec*square + L2_prec
				end if

			end do
		end do
		


		L2 = sqrt(L2)
		read(14, *), L_inf

		L2_prec = sqrt(L2_prec)
		L_inf_prec = MAXVAL(abs(this.surface_precise(-this.lon_max+1:this.lon_max-1, -this.lat_max+1:this.lat_max-1)))

			write(11, *),time*this.convert_time,"	", abs(L1/L1_prec)
			write(12, *),time*this.convert_time,"	", abs(L2/L2_prec)
			write(13, *),time*this.convert_time,"	", abs(L_inf/abs(L_inf_prec) - 1)


	end subroutine



end module