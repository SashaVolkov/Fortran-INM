module printer_ncdf

	use func_var, Only: f_var
	use grid_var, Only: g_var
	use diagnostic_mod, Only: diagnostic
	use netcdf

	implicit none
	include"mpif.h"

	Private
	Public :: printer

	Type printer

		Integer(4) :: Wid, Courantid

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: to_print => to_print
		Procedure, Public :: print_grid => print_grid
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS



	Subroutine init(this, dim, step, Tmax, speedup, time, grid_id, ncid, ncid_gr, rescale, grid_type)

		Class(printer) :: this
		Integer(4), intent(in) :: dim, step, Tmax, speedup, rescale, grid_type
		Integer(4), intent(out) :: time, grid_id, ncid, ncid_gr

		Integer(4) status, face, xid, yid, faceid, llid, gr_xid, gr_yid, gr_faceid, Wid, Courantid
		character(40) istring, istring1
		character(80) path1, path2

		write(istring, *) 2*dim
		write(istring1, *) 2*step

		istring = trim(adjustl(istring))//'/'//trim(adjustl(istring1))//'th'
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

		path1 = trim('datFiles/'//trim(adjustl(istring))//'surface.nc')
		path2 = trim('datFiles/'//trim(adjustl(istring))//'grid.nc')


		status = nf90_create (path = path1, cmode = IOR(NF90_NETCDF4,IOR(NF90_MPIIO,NF90_CLOBBER)),&
		 comm = MPI_COMM_WORLD, info = MPI_INFO_NULL, ncid = ncid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		status = nf90_def_dim (ncid, "x", 2*dim, xid)
		status = nf90_def_dim (ncid, "y", 2*dim, yid)
		status = nf90_def_dim (ncid, "face", 6, faceid)
		status = nf90_def_dim (ncid, "time", Tmax/speedup + 1, time)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		status = nf90_def_var (ncid, "Level", NF90_FLOAT, (/ xid, yid, faceid, time/), Wid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)
		status = nf90_enddef (ncid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		status = nf90_def_var (ncid, "CFL", NF90_FLOAT, (/ xid, yid, faceid, time/), Courantid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)
		status = nf90_enddef (ncid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)


		status = nf90_create (path = path2, cmode = IOR(NF90_NETCDF4,IOR(NF90_MPIIO,NF90_CLOBBER)),&
		 comm = MPI_COMM_WORLD, info = MPI_INFO_NULL, ncid = ncid_gr)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		status = nf90_def_dim (ncid_gr, "ll", 2, llid)
		status = nf90_def_dim (ncid_gr, "x", 2*dim, gr_xid)
		status = nf90_def_dim (ncid_gr, "y", 2*dim, gr_yid)
		status = nf90_def_dim (ncid_gr, "face", 6, gr_faceid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		status = nf90_def_var (ncid_gr, "latlon", NF90_DOUBLE, (/ llid, gr_xid, gr_yid, gr_faceid/), grid_id)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)
		status = nf90_enddef (ncid_gr)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		this.Courantid = Courantid;  this.Wid = Wid

	end Subroutine



	Subroutine to_print(this, var, diagn, time, speedup, ncid, id)

		Class(printer) :: this
		Class(f_var) :: var
		Class(diagnostic) :: diagn
		Integer(4), intent(in) :: time, speedup, ncid, id

		Integer(4) x, y, face, ier
		Integer(4) status, t, ns_y, ns_x, nf_y, nf_x, Ysize, Xsize, Wid, Courantid

		Courantid = this.Courantid;  Wid = this.Wid

		ns_y = var.ns_y;  nf_y = var.nf_y
		ns_x = var.ns_x;  nf_x = var.nf_x
		Ysize = 1 + nf_y - ns_y; Xsize = 1 + nf_x - ns_x

		t = 1+time/speedup

		do face = 1, 6
			status = nf90_put_var(ncid, Wid, real(var.h_height(ns_x:nf_x, ns_y:nf_y, face),4),&
			 start = (/ ns_x, ns_y, face, t/), count = (/ Xsize, Ysize, 1, 1/))
			if(status /= nf90_NoErr) print *, nf90_strerror(status) , id

			status = nf90_put_var(ncid, Courantid, real(diagn.CFL(ns_x:nf_x, ns_y:nf_y, face),4),&
			 start = (/ ns_x, ns_y, face, t/), count = (/ Xsize, Ysize, 1, 1/))
			if(status /= nf90_NoErr) print *, nf90_strerror(status) , id
		end do
	end Subroutine


	Subroutine print_grid(this, grid, grid_id, ncid_gr)

		Class(printer) :: this
		Class(g_var) :: grid
		Integer(4), intent(in) :: grid_id, ncid_gr
		Integer(4) x, y, face, ier, dim, status

		dim = grid.dim

			status = nf90_put_var(ncid_gr, grid_id, grid.latlon_c(1:2, 1:2*dim, 1:2*dim, 1:6),&
			 start = (/1, 1, 1, 1/), count = (/2, 2*dim, 2*dim, 6/))
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
	end Subroutine



	Subroutine deinit(this, ncid, ncid_gr)
		Class(printer) :: this
		Integer(4), intent(in) :: ncid, ncid_gr
		Integer(4) status

		status = nf90_close (ncid)
		status = nf90_close (ncid_gr)
	end Subroutine



end module