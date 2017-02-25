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

		Integer(4) :: Wid, Courantid, grid_id, ncid, ncid_gr, speedup, id

		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: to_print => to_print
		Procedure, Public :: print_grid => print_grid
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS



	Subroutine init(this, grid, Tmax, speedup, time, rescale, grid_type)

		Class(printer) :: this
		Class(g_var) :: grid
		Integer(4), intent(in) :: Tmax, speedup, rescale, grid_type
		Integer(4), intent(out) :: time

		Integer(4) status, face, xid, yid, dim, step, faceid, llid, gr_xid, gr_yid
		Integer(4) gr_faceid, Wid, Courantid, Precid, grid_id, ncid, ncid_gr, id, ier
		character(40) istring, istring1
		character(80) path1, path2, path3

		this.speedup = speedup;  dim = grid.dim;  step = grid.step - 1
		call MPI_Comm_rank(MPI_COMM_WORLD,id,ier)
		this.id = id

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
		path2 = trim('datFiles/'//trim(adjustl(istring))//'grid.dat')
		path3 = trim('datFiles/'//trim(adjustl(istring))//'square.dat')


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

		status = nf90_def_var (ncid, "Precise", NF90_FLOAT, (/ xid, yid, faceid, time/), Precid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)
		status = nf90_enddef (ncid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		this.Courantid = Courantid;  this.Wid = Wid; this.grid_id = grid_id; this.ncid = ncid; this.ncid_gr = ncid_gr

		if(id==0) then
			open(40,file=path2)
			open(42,file=path3)
		end if

	end Subroutine



	Subroutine to_print(this, var, diagn, time)

		Class(printer) :: this
		Class(f_var) :: var
		Class(diagnostic) :: diagn
		Integer(4), intent(in) :: time

		Integer(4) x, y, face, ier
		Integer(4) status, t, ns_y, ns_x, nf_y, nf_x, Ysize, Xsize, Wid, Courantid

		Courantid = this.Courantid;  Wid = this.Wid

		ns_y = var.ns_y;  nf_y = var.nf_y
		ns_x = var.ns_x;  nf_x = var.nf_x
		Ysize = 1 + nf_y - ns_y; Xsize = 1 + nf_x - ns_x

		t = 1+time/this.speedup

		do face = 1, 6
			status = nf90_put_var(this.ncid, Wid, real(var.h_height(ns_x:nf_x, ns_y:nf_y, face),4),&
			 start = (/ ns_x, ns_y, face, t/), count = (/ Xsize, Ysize, 1, 1/))
			if(status /= nf90_NoErr) print *, nf90_strerror(status) , this.id

			status = nf90_put_var(this.ncid, Courantid, real(diagn.CFL(ns_x:nf_x, ns_y:nf_y, face),4),&
			 start = (/ ns_x, ns_y, face, t/), count = (/ Xsize, Ysize, 1, 1/))
			if(status /= nf90_NoErr) print *, nf90_strerror(status) , this.id
		end do
	end Subroutine


	Subroutine print_grid(this, grid)

		Class(printer) :: this
		Class(g_var) :: grid
		Integer(4) x, y, face, ier, dim, status

		dim = grid.dim
		if(this.id==0) then
			write(40, *),grid.latlon_c(1:2,1:2*dim, 1:2*dim,1:6)
			write(42, *),real(grid.square(1:2*dim, 1:2*dim),4)
			close(40)
			close(42)
		end if
	end Subroutine


	Subroutine deinit(this)
		Class(printer) :: this
		Integer(4) :: status
		status = nf90_close (this.ncid)
	end Subroutine


end module