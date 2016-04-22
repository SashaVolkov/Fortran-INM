module printer_ncdf

	Use netcdf

	implicit none

	Private
	Public :: printer

	Type printer
		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: to_print => to_print
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS



	subroutine init(this)

		Class(printer) :: this
		integer(4) dim, Tmax, time, Wid, latitude, longitude, ncid

		integer(4) status
		character(40) name

		name = "one.nc"

			status = nf90_create (path = trim('/home/sasha/Fortran/Cubic_grid/solver/datFiles/'//name), &
	cmode = NF90_CLOBBER, ncid = ncid)

			status = nf90_def_dim (ncid, "latitude", dim, latitude)
			status = nf90_def_dim (ncid, "longitude", dim, longitude)
			status = nf90_def_dim (ncid, "time", Tmax, time)

			status = nf90_def_var (ncid, "water", NF90_REAL, (/ latitude, longitude, time/), Wid)
			status = nf90_enddef (ncid)

	end subroutine



	! 		call m.to_print(fprev, g, t, name, Tmax, Wid, xid, yid, time, ncid)


		subroutine to_print(this)

			Class(printer) :: this

	! 		Integer x, y, i, j
	! 		Integer status 
	! 		Integer(4), Intent(In) :: t, Tmax, Wid, xid, yid, time, ncid
	! 		character(40), Intent(In) :: name

	! 		Real(8) W_mass(g.ns_y:g.nf_y, g.ns_x:g.nf_x)

	! 		do j = g.ns_x, g.nf_x
	! 			do i = g.ns_y, g.nf_y
	! 				W_mass(i, j) = f.d(i, j)
	! 			end do
	! 		end do


	! 		! ! указатель на первый элемент массива varval , число элементов
	! 		  status = nf90_put_var (ncid, Wid, W_mass, (/ g.ns_y, g.ns_x, t/), (/ g.nf_y - g.ns_y + 1, g.nf_x - g.ns_x + 1, 1/) )

			end subroutine

	! function nf90_put_var(ncid, varid, values, start, count, stride, map)

	subroutine deinit(this)
		Class(printer) :: this
		integer(4) status, ncid

		status = nf90_close (ncid)
	end subroutine



end module