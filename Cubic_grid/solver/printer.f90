module printer_ncdf

	use func_var, Only: f_var
	use netcdf

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



	subroutine init(this, dim, Tmax, speedup, time, Wid, xid, yid, ncid)

		Class(printer) :: this
		integer(4), intent(in) :: dim, Tmax, speedup
		integer(4), intent(out) :: time, Wid, xid, yid, ncid(1:6)

		integer(4) status, face_index
		character(40) istring


		do face_index = 1, 6
		write(istring(1:1), '(i1.1)') face_index
		istring = istring//'.nc'

			status = nf90_create (path = trim('/home/sasha/Fortran/Cubic_grid/solver/datFiles/'//"face"//istring), &
	cmode = NF90_CLOBBER, ncid = ncid(face_index))
			if(status /= nf90_NoErr) print *, nf90_strerror(status)


			status = nf90_def_dim (ncid(face_index), "x", 2*dim+1, xid)
			status = nf90_def_dim (ncid(face_index), "y", 2*dim+1, yid)
			status = nf90_def_dim (ncid(face_index), "time", Tmax/speedup + 1, time)
			if(status /= nf90_NoErr) print *, nf90_strerror(status)

			status = nf90_def_var (ncid(face_index), "water", NF90_DOUBLE, (/ xid, yid, time/), Wid)
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
			status = nf90_enddef (ncid(face_index))
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
			end do

	end subroutine



	! 		call m.to_print(fprev, g, t, name, Tmax, Wid, xid, yid, time, ncid)


	subroutine to_print(this, var, dim, time, speedup, Wid, ncid)

		Class(printer) :: this
		Class(f_var) :: var
		integer(4), intent(in) :: dim, time, speedup, Wid, ncid(1:6)

		integer(4) i, j, face_index
		integer(4) status


		Real(8) W_mass( -dim:dim, -dim:dim)

		do face_index = 1, 6
			do j = -dim, dim
				do i = -dim, dim
					W_mass(i, j) = var.h_height(i, j, face_index)
				end do
			end do

			status = nf90_put_var(ncid(face_index), Wid, W_mass, start = (/ 1, 1, 1+time/speedup/), count = (/ 2*dim+1, 2*dim+1, 1/))
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
		end do

			! print '("  h = ", f10.2, " m")', W_mass(0, 0)

		end subroutine



	subroutine deinit(this)
		Class(printer) :: this
		integer(4) status, ncid

		status = nf90_close (ncid)
	end subroutine



end module