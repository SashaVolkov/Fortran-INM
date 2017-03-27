module scan_print

	use netcdf

	implicit none

	Private
	Public :: printer


	Type printer

	integer(4) :: dim, step, ncid, ncid_to, grid_id, Wid(8), Wid_to(7)
	integer(4) :: lon_max, lat_max, nc_or_dat, point_id, ncid_point, point_find
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



	subroutine init(this, dim, step, all_time, convert_time, rescale, grid_type, lon_max, lat_max, nc_or_dat)

		Class(printer) :: this
		integer(4), intent(in) :: dim, step, all_time, rescale, grid_type, lon_max, lat_max, nc_or_dat

		integer(4) status, face, ncid, ncid_to, nvars, grid_id(1:1), Wid(8), Wid_to(7), time, Courantid_to, point_id(1:1), ncid_point, point_find
		logical file_exist
		integer(4) latid, lonid, coord
		real(8) convert_time
		character(80) istring, istring1
		character(80) path1, path2, path3, path4

		this.dim = dim;  this.lon_max = lon_max;  this.lat_max = lat_max;  this.convert_time = convert_time
		this.nc_or_dat = nc_or_dat;  this.step = 2

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
		path2 = trim(trim(adjustl(istring))//'grid.dat')
		if(this.nc_or_dat == 0) then
			path3 = trim(trim(adjustl(istring))//'surface_ll.nc')
		else
			path3 = trim(trim(adjustl(istring))//'surface_ll.dat')
		end if
		path4 = trim(trim(adjustl(istring))//'closest_point_ll_C.dat')

		open(40, file=path2)


		status = nf90_open(path = path1, mode = NF90_WRITE, ncid = ncid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status), "init1", path1
		status = nf90_inq_varids(ncid, nvars, Wid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status), "init2", path1

		if(this.nc_or_dat == 0) then
			status = nf90_create(path = path3, cmode = NF90_CLOBBER, ncid = ncid_to)
			status = nf90_def_dim (ncid_to, "lon", 2*this.lon_max+1, lonid)
			status = nf90_def_dim (ncid_to, "lat", 2*this.lat_max+1, latid)
			status = nf90_def_dim (ncid_to, "time", all_time, time)
			status = nf90_def_var (ncid_to, "Level", NF90_FLOAT, (/ lonid, latid, time/), Wid_to(1))
			status = nf90_def_var (ncid_to, "Lon_vel", NF90_FLOAT, (/ lonid, latid, time/), Wid_to(2))
			status = nf90_def_var (ncid_to, "Lat_vel", NF90_FLOAT, (/ lonid, latid, time/), Wid_to(3))
			status = nf90_def_var (ncid_to, "Level_er", NF90_FLOAT, (/ lonid, latid, time/), Wid_to(4))
			status = nf90_def_var (ncid_to, "Lon_vel_er", NF90_FLOAT, (/ lonid, latid, time/), Wid_to(5))
			status = nf90_def_var (ncid_to, "Lat_vel_er", NF90_FLOAT, (/ lonid, latid, time/), Wid_to(6))
			status = nf90_def_var (ncid_to, "Courant", NF90_FLOAT, (/ lonid, latid, time/), Wid_to(7))
			status = nf90_enddef (ncid_to)
			if(status /= nf90_NoErr) print *, nf90_strerror(status), "init3"
			this.ncid_to = ncid_to;  this.Wid_to = Wid_to
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

		this.ncid = ncid;  this.Wid = Wid;  this.point_find = point_find

	end subroutine



	subroutine scan_surf(this, time, surface_off)

		Class(printer) :: this
		integer(4), intent(in) :: time
		real(8), intent(out) :: surface_off(1-this.step:2*this.dim+this.step, 1-this.step:2*this.dim+this.step, 6, 7)
		integer(4) x, y, face, ier, status, ncid, Wid(8), dim

		dim = this.dim;  ncid = this.ncid;  Wid = this.Wid

		status = nf90_get_var(ncid, Wid(1), surface_off(1:2*dim, 1:2*dim, 1:6, 1),&
		 start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))

		status = nf90_get_var(ncid, Wid(2), surface_off(1:2*dim, 1:2*dim, 1:6, 2),&
		 start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))

		status = nf90_get_var(ncid, Wid(3), surface_off(1:2*dim, 1:2*dim, 1:6, 3),&
		 start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))

		status = nf90_get_var(ncid, Wid(4), surface_off(1:2*dim, 1:2*dim, 1:6, 4),&
		 start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))

		status = nf90_get_var(ncid, Wid(5), surface_off(1:2*dim, 1:2*dim, 1:6, 5),&
		 start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))

		status = nf90_get_var(ncid, Wid(6), surface_off(1:2*dim, 1:2*dim, 1:6, 6),&
		 start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))

		status = nf90_get_var(ncid, Wid(7), surface_off(1:2*dim, 1:2*dim, 1:6, 7),&
		 start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))
		if(status /= nf90_NoErr) print *, nf90_strerror(status), "scan_surf"
	end subroutine




	subroutine print_surf_prec_cube(this, time, precise, surface_to)

		Class(printer) :: this
		integer(4), intent(in) :: time
		real(4), intent(in) :: precise(1:2*this.dim, 1:2*this.dim, 1:6)
		real(8), intent(in) :: surface_to(1-this.step:2*this.dim+this.step, 1-this.step:2*this.dim+this.step, 7)
		integer(4) status, ncid, dim

		dim = this.dim;  ncid = this.ncid

		! status = nf90_put_var(ncid, this.Wid(8), precise(1:2*dim, 1:2*dim, 1:6)-real(surface_to(1:2*dim, 1:2*dim, 1:6),4),&
		!  start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))
		! if(status /= nf90_NoErr) print *, nf90_strerror(status), "print_surf_cube"
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
		real(8), intent(out) :: grid(2, 1-this.step:2*this.dim+this.step, 1-this.step:2*this.dim+this.step, 7)

		read(40, *), grid(1:2, 1:2*this.dim, 1:2*this.dim, 1:6)
		close(40)

	end subroutine



	subroutine scan_point(this, point)
		Class(printer) :: this
		integer(4), intent(out) :: point(1:3, -this.lon_max:this.lon_max, -this.lat_max:this.lat_max)
		read(16, *) point(1:3, -this.lon_max:this.lon_max, -this.lat_max:this.lat_max)
	end subroutine



	subroutine print_surf(this, surface_to, time)
		Class(printer) :: this
		real(8), intent(in) :: surface_to(-this.lon_max:this.lon_max, -this.lat_max:this.lat_max, 7)
		! real(4), intent(in) :: surface_precise(-this.lon_max:this.lon_max, -this.lat_max:this.lat_max)
		integer(4), intent(in) :: time
		integer(4) status, Wid_to(7), ncid_to, Courantid_to

		ncid_to = this.ncid_to;  Wid_to = this.Wid_to

		if(this.nc_or_dat == 0) then
			status = nf90_put_var(ncid_to, Wid_to(1), real(surface_to(:,:,1),4),&
			 start = (/1, 1, time/), count = (/2*this.lon_max+1, 2*this.lat_max+1, 1/))

			status = nf90_put_var(ncid_to, Wid_to(2), real(surface_to(:,:,2),4),&
			 start = (/1, 1, time/), count = (/2*this.lon_max+1, 2*this.lat_max+1, 1/))

			status = nf90_put_var(ncid_to, Wid_to(3), real(surface_to(:,:,3),4),&
			 start = (/1, 1, time/), count = (/2*this.lon_max+1, 2*this.lat_max+1, 1/))

			status = nf90_put_var(ncid_to, Wid_to(4), real(surface_to(:,:,4),4),&
			 start = (/1, 1, time/), count = (/2*this.lon_max+1, 2*this.lat_max+1, 1/))

			status = nf90_put_var(ncid_to, Wid_to(5), real(surface_to(:,:,5),4),&
			 start = (/1, 1, time/), count = (/2*this.lon_max+1, 2*this.lat_max+1, 1/))

			status = nf90_put_var(ncid_to, Wid_to(6), real(surface_to(:,:,6),4),&
			 start = (/1, 1, time/), count = (/2*this.lon_max+1, 2*this.lat_max+1, 1/))

			status = nf90_put_var(ncid_to, Wid_to(7), real(surface_to(-this.lon_max:this.lon_max, -this.lat_max:this.lat_max, 7),4),&
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
		integer(4) :: ncid, ncid_to, status
		ncid = this.ncid; ncid_to = this.ncid_to

		status = nf90_close(ncid)
		if(this.nc_or_dat == 0) then
			status = nf90_close(ncid_to)
		else
			close(15)
		end if
		close(40)
	end subroutine



end module