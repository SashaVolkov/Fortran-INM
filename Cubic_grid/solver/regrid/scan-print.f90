module scan_print

	use netcdf

	implicit none

	Private
	Public :: printer


	Type printer

	integer(4) :: dim, step, ncid, ncid_gr, ncid_to, grid_id, Wid, Wid_to
	integer(4) :: lon_max, lat_max, nc_or_dat, Courantid, Precid, Courantid_to, point_id, ncid_point, point_find
	real(8) :: convert_time
		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: scan_surf => scan_surf
		Procedure, Public :: scan_precise => scan_precise
		Procedure, Public :: scan_grid => scan_grid
		Procedure, Public :: scan_point => scan_point
		Procedure, Public :: print_surf => print_surf
		Procedure, Public :: print_surf_prec_cube => print_surf_prec_cube
		Procedure, Public :: print_point => print_point
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS



	subroutine init(this, dim, step, all_time, convert_time, rescale, grid_type, lon_max, lat_max)

		Class(printer) :: this
		integer(4), intent(in) :: dim, step, all_time, rescale, grid_type, lon_max, lat_max

		integer(4) status, face, ncid, ncid_to, ncid_gr, nvars, grid_id(1:1), Wid(1:3), Wid_to, time, Courantid_to, point_id(1:1), ncid_point, point_find
		logical file_exist
		integer(4) latid, lonid, coord
		real(8) convert_time
		character(80) istring, istring1
		character(80) path1, path2, path3, path4

		this.dim = dim;  this.lon_max = lon_max;  this.lat_max = lat_max;  this.convert_time = convert_time
		this.nc_or_dat = 0;  this.step = 2

		write(istring, *) 2*dim
		write(istring1, *) 2*step

		istring = '../datFiles/'//trim(adjustl(istring))//'/'//trim(adjustl(istring1))//'th'

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

		path1 = trim(trim(adjustl(istring))//'surface.nc')
		path2 = trim(trim(adjustl(istring))//'grid.nc')
		if(this.nc_or_dat == 0) then
			path3 = trim(trim(adjustl(istring))//'surface_ll.nc')
		else
			path3 = trim(trim(adjustl(istring))//'surface_ll.dat')
		end if
		path4 = trim(trim(adjustl(istring))//'closest_point_ll_C.dat')



		status = nf90_open(path = path1, mode = NF90_WRITE, ncid = ncid)
		status = nf90_inq_varids(ncid, nvars, Wid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		status = nf90_open (path = path2, mode = NF90_NOWRITE, ncid = ncid_gr)
		status = nf90_inq_varids(ncid_gr, nvars, grid_id)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		if(this.nc_or_dat == 0) then
			status = nf90_create(path = path3, cmode = NF90_CLOBBER, ncid = ncid_to)
			status = nf90_def_dim (ncid_to, "lon", 2*this.lon_max+1, lonid)
			status = nf90_def_dim (ncid_to, "lat", 2*this.lat_max+1, latid)
			status = nf90_def_dim (ncid_to, "time", all_time, time)
			status = nf90_def_var (ncid_to, "Level", NF90_FLOAT, (/ lonid, latid, time/), Wid_to)
			status = nf90_def_var (ncid_to, "Courant", NF90_FLOAT, (/ lonid, latid, time/), Courantid_to)
			status = nf90_enddef (ncid_to)
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
			this.ncid_to = ncid_to;  this.Wid_to = Wid_to;  this.Courantid_to = Courantid_to
		else
			open(15,file=path3,access="direct",recl=(2*this.lat_max+1)*(2*this.lon_max+1))
		end if

			inquire(file=path4,exist=file_exist)
			if(file_exist) then
				point_find = 0
				open(UNIT=16,file=path4,FORM="FORMATTED",STATUS="OLD",ACTION="READ")
			else
				point_find = 1
				open(UNIT=16,file=path4,FORM="FORMATTED",STATUS="NEW",ACTION="WRITE")
			end if

		this.ncid = ncid;  this.ncid_gr = ncid_gr;  this.grid_id = grid_id(1); this.Wid = Wid(1); this.Courantid = Wid(2);  this.Precid = Wid(3);  this.point_find = point_find
		! this.ncid_point = ncid_point;  this.point_id = point_id(1)

	end subroutine



	subroutine scan_surf(this, time, surface_off)

		Class(printer) :: this
		integer(4), intent(in) :: time
		real(8), intent(out) :: surface_off(1-this.step:2*this.dim+this.step, 1-this.step:2*this.dim+this.step, 1:6, 1:2)
		integer(4) x, y, face, ier, status, ncid, Wid, dim, Courantid

		dim = this.dim;  ncid = this.ncid;  Wid = this.Wid;  Courantid = this.Courantid

		status = nf90_get_var(ncid, Wid, surface_off(1:2*dim, 1:2*dim, 1:6, 1),&
		 start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))

		status = nf90_get_var(ncid, Courantid, surface_off(1:2*dim, 1:2*dim, 1:6, 2),&
		 start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))
		if(status /= nf90_NoErr) print *, nf90_strerror(status), "scan_surf"
	end subroutine




	subroutine print_surf_prec_cube(this, time, surface_to)

		Class(printer) :: this
		integer(4), intent(in) :: time
		real(4), intent(in) :: surface_to(1:2*this.dim, 1:2*this.dim, 1:6)
		integer(4) status, ncid, dim, Precid

		dim = this.dim;  ncid = this.ncid;  Precid = this.Precid

		status = nf90_put_var(ncid, Precid, surface_to(1:2*dim, 1:2*dim, 1:6),&
		 start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))
		if(status /= nf90_NoErr) print *, nf90_strerror(status), "print_surf_cube"
	end subroutine




	subroutine scan_precise(this, time, surface_precise)

		Class(printer) :: this
		integer(4), intent(in) :: time
		real(4), intent(out) :: surface_precise(-this.lon_max:this.lon_max, -this.lat_max:this.lat_max)

		open(22,file="../datFiles/f.dat",access="direct", recl=(2*this.lon_max+1)*(2*this.lat_max+1))
		read(22, rec=time) surface_precise(:,:)
		close(22)
	end subroutine



	subroutine scan_grid(this, grid)

		Class(printer) :: this
		real(8), intent(out) :: grid(1:2, 1-this.step:2*this.dim+this.step, 1-this.step:2*this.dim+this.step, 1:6)
		integer(4) x, y, face, ier, dim, status, grid_id, ncid_gr

		dim = this.dim;  ncid_gr = this.ncid_gr;  grid_id = this.grid_id

		status = nf90_get_var(ncid_gr, grid_id, grid(1:2, 1:2*dim, 1:2*dim, 1:6),&
		 start = (/1, 1, 1, 1/), count = (/2, 2*dim, 2*dim, 6/))
		if(status /= nf90_NoErr) print *, nf90_strerror(status), "scan_grid"
	end subroutine



	subroutine scan_point(this, point)
		Class(printer) :: this
		integer(4), intent(out) :: point(1:3, -this.lon_max:this.lon_max, -this.lat_max:this.lat_max)
		read(16, *) point(1:3, -this.lon_max:this.lon_max, -this.lat_max:this.lat_max)
	end subroutine



	subroutine print_surf(this, surface_to, surface_precise, time)
		Class(printer) :: this
		real(8), intent(in) :: surface_to(-this.lon_max:this.lon_max, -this.lat_max:this.lat_max, 2)
		real(4), intent(in) :: surface_precise(-this.lon_max:this.lon_max, -this.lat_max:this.lat_max)
		integer(4), intent(in) :: time
		integer(4) status, Wid_to, ncid_to, Courantid_to

		ncid_to = this.ncid_to;  Wid_to = this.Wid_to;  Courantid_to = this.Courantid_to

		if(this.nc_or_dat == 0) then
			status = nf90_put_var(ncid_to, Wid_to, real(surface_to(:,:,1),4),&
			 start = (/1, 1, time/), count = (/2*this.lon_max+1, 2*this.lat_max+1, 1/))

			status = nf90_put_var(ncid_to, Courantid_to, real(surface_to(-this.lon_max:this.lon_max, -this.lat_max:this.lat_max, 2),4),&
			 start = (/1, 1, time/), count = (/2*this.lon_max+1, 2*this.lat_max+1, 1/))
			if(status /= nf90_NoErr) print *, nf90_strerror(status)
		else
			write(15, rec=time) real(surface_to(:,:, 1),4)
		end if
	end subroutine


	subroutine print_point(this, point)
		Class(printer) :: this
		integer(4), intent(in) :: point(1:3, -this.lon_max:this.lon_max, -this.lat_max:this.lat_max)
			write(16, *) point(1:3, -this.lon_max:this.lon_max, -this.lat_max:this.lat_max)
	end subroutine



	subroutine deinit(this)
		Class(printer) :: this
		integer(4) :: ncid, ncid_to, ncid_gr, status
		ncid_gr = this.ncid_gr;  ncid = this.ncid; ncid_to = this.ncid_to

		status = nf90_close(ncid)
		status = nf90_close(ncid_gr)
		if(this.nc_or_dat == 0) then
			status = nf90_close(ncid_to)
		else
			close(15)
		end if
	end subroutine



end module