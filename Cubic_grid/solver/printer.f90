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

		integer(4) :: Wid, Courantid

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: to_print => to_print
		Procedure, Public :: print_grid => print_grid
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS



	subroutine init(this, dim, Tmax, speedup, time, grid_id, ncid, ncid_gr, rescale, grid_type)

		Class(printer) :: this
		integer(4), intent(in) :: dim, Tmax, speedup, rescale, grid_type
		integer(4), intent(out) :: time, grid_id, ncid, ncid_gr

		integer(4) status, face, xid, yid, faceid, llid, gr_xid, gr_yid, gr_faceid, Wid, Courantid
		character(40) istring
		character(80) path1, path2

		write(istring, *) 2*dim

		if(grid_type == 0) then
			if(rescale == 0) then
				path1 = trim('datFiles/'//trim(adjustl(istring))//'/simple/surface.nc')
				path2 = trim('datFiles/'//trim(adjustl(istring))//'/simple/grid.nc')
			else if(rescale == 1) then
				path1 = trim('datFiles/'//trim(adjustl(istring))//'/tan/surface.nc')
				path2 = trim('datFiles/'//trim(adjustl(istring))//'/tan/grid.nc')
			else if(rescale == 2) then
				path1 = trim('datFiles/'//trim(adjustl(istring))//'/exp/surface.nc')
				path2 = trim('datFiles/'//trim(adjustl(istring))//'/exp/grid.nc')
			end if
		else if(grid_type == 1) then
			path1 = trim('datFiles/'//trim(adjustl(istring))//'/equiang/surface.nc')
			path2 = trim('datFiles/'//trim(adjustl(istring))//'/equiang/grid.nc')
		end if

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

	end subroutine



	subroutine to_print(this, var, diagn, time, speedup, ncid, id)

		Class(printer) :: this
		Class(f_var) :: var
		Class(diagnostic) :: diagn
		integer(4), intent(in) :: time, speedup, ncid, id

		integer(4) x, y, face, ier
		integer(4) status, t, ns_y, ns_x, nf_y, nf_x, Ysize, Xsize, Wid, Courantid

		Courantid = this.Courantid;  Wid = this.Wid

		ns_y = var.ns_y;  nf_y = var.nf_y
		ns_x = var.ns_x;  nf_x = var.nf_x
		Ysize = var.Ysize; Xsize = var.Xsize

		t = 1+time/speedup

		do face = 1, 6
			status = nf90_put_var(ncid, Wid, real(var.h_height(ns_x:nf_x, ns_y:nf_y, face),4),&
			 start = (/ ns_x, ns_y, face, t/), count = (/ Xsize, Ysize, 1, 1/))
			if(status /= nf90_NoErr) print *, nf90_strerror(status) , id

			status = nf90_put_var(ncid, Courantid, real(diagn.CFL(ns_x:nf_x, ns_y:nf_y, face),4),&
			 start = (/ ns_x, ns_y, face, t/), count = (/ Xsize, Ysize, 1, 1/))
			if(status /= nf90_NoErr) print *, nf90_strerror(status) , id
		end do
	end subroutine


	subroutine print_grid(this, grid, grid_id, ncid_gr)

		Class(printer) :: this
		Class(g_var) :: grid
		integer(4), intent(in) :: grid_id, ncid_gr
		integer(4) x, y, face, ier, dim, status

		dim = grid.dim

			status = nf90_put_var(ncid_gr, grid_id, grid.latlon_c(1:2, 1:2*dim, 1:2*dim, 1:6),&
			 start = (/1, 1, 1, 1/), count = (/2, 2*dim, 2*dim, 6/))
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
	end subroutine



	subroutine deinit(this, ncid, ncid_gr)
		Class(printer) :: this
		integer(4), intent(in) :: ncid, ncid_gr
		integer(4) status

		status = nf90_close (ncid)
		status = nf90_close (ncid_gr)
	end subroutine



end module