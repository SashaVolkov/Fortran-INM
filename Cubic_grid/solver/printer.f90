module printer_ncdf


	implicit none

	Private
	Public :: printer

	Type printer
		CONTAINS
		Procedure, Public :: init => init
	End Type

	CONTAINS



	subroutine init(this)

		Class(printer) :: this

	! 		status = nf90_create (path = trim('/home/sasha/Fortran/Shallow_Water/datFiles/'//name), &
	! cmode = IOR(NF90_NETCDF4,IOR(NF90_MPIIO,NF90_CLOBBER)), comm = MPI_COMM_WORLD, info = MPI_INFO_NULL, ncid = ncid)

	! 		status = nf90_def_dim (ncid, "x", g.StepsX, xid)
	! 		status = nf90_def_dim (ncid, "y", g.StepsY, yid)
	! 		status = nf90_def_dim (ncid, "t", Tmax, tid)

	! 		status = nf90_def_var (ncid, "water", NF90_REAL, (/ yid, xid, tid/), Wid)
	! 		status = nf90_enddef  (ncid)
	! 		call m.to_print(fprev, g, t, name, Tmax, Wid, xid, yid, tid, ncid)


	end subroutine



end module