module printer_ncdf

	use func_var, Only: f_var
	use netcdf

	implicit none
	include"mpif.h"

	Private
	Public :: printer

	Type printer
		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: to_print => to_print
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS



	subroutine init(this, dim, Tmax, speedup, time, Wid, xid, yid, ncid, rescale)

		Class(printer) :: this
		integer(4), intent(in) :: dim, Tmax, speedup, rescale
		integer(4), intent(out) :: time, Wid, xid, yid, ncid(1:6)

		integer(4) status, face
		character(40) istring


		do face = 1, 6
		write(istring(1:1), '(i1.1)') face
		istring = istring//'.nc'

		if(rescale == 0) then
			status = nf90_create (path = trim('datFiles/simple/'//"face"//istring), &
	cmode = IOR(NF90_NETCDF4,IOR(NF90_MPIIO,NF90_CLOBBER)), comm = MPI_COMM_WORLD, info = MPI_INFO_NULL, ncid = ncid(face))
		else if(rescale == 1) then
			status = nf90_create (path = trim('datFiles/tan/'//"face"//istring), &
	cmode = IOR(NF90_NETCDF4,IOR(NF90_MPIIO,NF90_CLOBBER)), comm = MPI_COMM_WORLD, info = MPI_INFO_NULL, ncid = ncid(face))
		else if(rescale == 2) then
			status = nf90_create (path = trim('datFiles/exp/'//"face"//istring), &
	cmode = IOR(NF90_NETCDF4,IOR(NF90_MPIIO,NF90_CLOBBER)), comm = MPI_COMM_WORLD, info = MPI_INFO_NULL, ncid = ncid(face))
		end if

			if(status /= nf90_NoErr) print *, nf90_strerror(status)


			status = nf90_def_dim (ncid(face), "x", 2*dim+1, xid)
			status = nf90_def_dim (ncid(face), "y", 2*dim+1, yid)
			status = nf90_def_dim (ncid(face), "time", Tmax/speedup + 1, time)
			if(status /= nf90_NoErr) print *, nf90_strerror(status)

			status = nf90_def_var (ncid(face), "water", NF90_DOUBLE, (/ xid, yid, time/), Wid)
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
			status = nf90_enddef (ncid(face))
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
			end do

	end subroutine



	! 		call m.to_print(fprev, g, t, name, Tmax, Wid, xid, yid, time, ncid)


	subroutine to_print(this, var, dim, time, speedup, Wid, ncid)

		Class(printer) :: this
		Class(f_var) :: var(1:6)
		integer(4), intent(in) :: dim, time, speedup, Wid, ncid(1:6)

		integer(4) i, j, face
		integer(4) status


		Real(8) W_mass( -dim:dim, -dim:dim)

		do face = 1, 6
			do j = -dim, dim
				do i = -dim, dim
					W_mass(i, j) = var(face).h_height(i, j)
				end do
			end do

			status = nf90_put_var(ncid(face), Wid, W_mass, start = (/ 1, 1, 1+time/speedup/), count = (/ 2*dim+1, 2*dim+1, 1/))
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