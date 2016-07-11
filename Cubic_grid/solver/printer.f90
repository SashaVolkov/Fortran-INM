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



	subroutine init(this, dim, Tmax, speedup, time, Wid, xid, yid, faceid, ncid, rescale)

		Class(printer) :: this
		integer(4), intent(in) :: dim, Tmax, speedup, rescale
		integer(4), intent(out) :: time, Wid, xid, yid, faceid, ncid

		integer(4) status, face
		character(40) istring


! 		do face = 1, 6
! 		write(istring(1:1), '(i1.1)') face
! 		istring = istring

		if(rescale == 0) then
			status = nf90_create (path = trim('datFiles/simple/'//"surface.nc"), &
	cmode = IOR(NF90_NETCDF4,IOR(NF90_MPIIO,NF90_CLOBBER)), comm = MPI_COMM_WORLD, info = MPI_INFO_NULL, ncid = ncid)
		else if(rescale == 1) then
			status = nf90_create (path = trim('datFiles/tan/'//"surface.nc"), &
	cmode = IOR(NF90_NETCDF4,IOR(NF90_MPIIO,NF90_CLOBBER)), comm = MPI_COMM_WORLD, info = MPI_INFO_NULL, ncid = ncid)
		else if(rescale == 2) then
			status = nf90_create (path = trim('datFiles/exp/'//"surface.nc"), &
	cmode = IOR(NF90_NETCDF4,IOR(NF90_MPIIO,NF90_CLOBBER)), comm = MPI_COMM_WORLD, info = MPI_INFO_NULL, ncid = ncid)
		end if

			if(status /= nf90_NoErr) print *, nf90_strerror(status)


			status = nf90_def_dim (ncid, "x", 2*dim, xid)
			status = nf90_def_dim (ncid, "y", 2*dim, yid)
			status = nf90_def_dim (ncid, "face", 6, faceid)
			status = nf90_def_dim (ncid, "time", Tmax/speedup + 1, time)
			if(status /= nf90_NoErr) print *, nf90_strerror(status)

			status = nf90_def_var (ncid, "water", NF90_DOUBLE, (/ xid, yid, faceid, time/), Wid)
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
			status = nf90_enddef (ncid)
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
! 			end do

	end subroutine



	! 		call m.to_print(fprev, g, t, name, Tmax, Wid, xid, yid, time, ncid)


	subroutine to_print(this, var, time, speedup, Wid, ncid, id)

		Class(printer) :: this
		Class(f_var) :: var
		integer(4), intent(in) :: time, speedup, Wid, ncid, id

		integer(4) x, y, face, ier
		integer(4) status, t, ns_y, ns_x, nf_y, nf_x, Ysize, Xsize
		real(8) W_mass(var.ns_x:var.nf_x, var.ns_y:var.nf_y)

		ns_y = var.ns_y;  nf_y = var.nf_y
		ns_x = var.ns_x;  nf_x = var.nf_x
		Ysize = var.Ysize; Xsize = var.Xsize

		t = 1+time/speedup


		do face = 1, 6

			status = nf90_put_var(ncid, Wid, var.h_height(ns_x:nf_x, ns_y:nf_y, face),&
			 start = (/ ns_x, ns_y, face, t/), count = (/ Xsize, Ysize, 1, 1/))

			if(status /= nf90_NoErr) print *, nf90_strerror(status) , id
		end do
	end subroutine



	subroutine deinit(this)
		Class(printer) :: this
		integer(4) status, ncid

		status = nf90_close (ncid)
	end subroutine



end module