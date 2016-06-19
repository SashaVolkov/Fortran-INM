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


			status = nf90_def_dim (ncid(face), "x", 2*dim, xid)
			status = nf90_def_dim (ncid(face), "y", 2*dim, yid)
			status = nf90_def_dim (ncid(face), "time", Tmax/speedup + 1, time)
			if(status /= nf90_NoErr) print *, nf90_strerror(status)

			status = nf90_def_var (ncid(face), "water", NF90_DOUBLE, (/ xid, yid, time/), Wid)
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
			status = nf90_enddef (ncid(face))
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
			end do

	end subroutine



	! 		call m.to_print(fprev, g, t, name, Tmax, Wid, xid, yid, time, ncid)


	subroutine to_print(this, var, time, speedup, Wid, ncid, id)

		Class(printer) :: this
		Class(f_var) :: var
		integer(4), intent(in) :: time, speedup, Wid, ncid(1:6), id

		integer(4) x, y, face, ier
		integer(4) status, t, ns_y, ns_x, nf_y, nf_x, Ysize, Xsize
		real(8) W_mass(var.ns_x:var.nf_x, var.ns_y:var.nf_y)

		ns_y = var.ns_y;  nf_y = var.nf_y
		ns_x = var.ns_x;  nf_x = var.nf_x
		Ysize = var.Ysize; Xsize = var.Xsize

		t = 1+time/speedup


		do face = 1, 6

			status = nf90_put_var(ncid(face), Wid, var.h_height(ns_x:nf_x, ns_y:nf_y, face),&
			 start = (/ ns_x, ns_y, t/), count = (/ Xsize, Ysize, 1/))

			if(status /= nf90_NoErr) print *, nf90_strerror(status) , id
		end do

	end subroutine



	subroutine deinit(this)
		Class(printer) :: this
		integer(4) status, ncid

		status = nf90_close (ncid)
	end subroutine



end module