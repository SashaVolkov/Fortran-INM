module scan_print

	use netcdf

	implicit none

	Private
	Public :: printer


	Type printer

	integer(4) :: dim, ncid, ncid_gr, ncid_to, grid_id, Wid, Wid_to, lon_max, lat_max
		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: scan_surf => scan_surf
		Procedure, Public :: scan_grid => scan_grid
		Procedure, Public :: print_surf => print_surf
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS



	subroutine init(this, dim, all_time, rescale, grid_type)

		Class(printer) :: this
		integer(4), intent(in) :: dim, all_time, rescale, grid_type

		integer(4) status, face, ncid, ncid_to, ncid_gr, nvars, grid_id(1:1), Wid(1:1), Wid_to, time
		integer(4) latid, lonid
		character(40) istring
		character(80) path1, path2, path3

		this.dim = dim;  this.lon_max = 180;  this.lat_max = 90

		write(istring, *) 2*dim

		if(grid_type == 0) then
			if(rescale == 0) then
				path1 = trim('../datFiles/simple/'//"surface_C"//trim(adjustl(istring))//".nc")
				path2 = trim('../datFiles/simple/'//"grid_C"//trim(adjustl(istring))//".nc")
				path3 = trim('../datFiles/simple/'//"surface_ll_C"//trim(adjustl(istring))//".dat")
			else if(rescale == 1) then
				path1 = trim('../datFiles/tan/'//"surface_C"//trim(adjustl(istring))//".nc")
				path2 = trim('../datFiles/tan/'//"grid_C"//trim(adjustl(istring))//".nc")
				path3 = trim('../datFiles/tan/'//"surface_ll_C"//trim(adjustl(istring))//".dat")
			else if(rescale == 2) then
				path1 = trim('../datFiles/'//"surface_conf_exp_C"//trim(adjustl(istring))//".nc")
				path2 = trim('../datFiles/'//"grid_conf_exp_C"//trim(adjustl(istring))//".nc")
				path3 = trim('../datFiles/'//"surface_conf_exp_ll_C"//trim(adjustl(istring))//".dat")
			end if
		else if(grid_type == 1) then
				path1 = trim('../datFiles/equiang/'//"surface_C"//trim(adjustl(istring))//".nc")
				path2 = trim('../datFiles/equiang/'//"grid_C"//trim(adjustl(istring))//".nc")
				path3 = trim('../datFiles/equiang/'//"surface_ll_C"//trim(adjustl(istring))//".dat")
		end if


		status = nf90_open(path = path1, mode = NF90_NOWRITE, ncid = ncid)
		status = nf90_inq_varids(ncid, nvars, Wid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		status = nf90_open (path = path2, mode = NF90_NOWRITE, ncid = ncid_gr)
		status = nf90_inq_varids(ncid_gr, nvars, grid_id)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)


		! status = nf90_create(path = path3, cmode = NF90_CLOBBER, ncid = ncid_to)
		! status = nf90_def_dim (ncid_to, "lon", 2*this.lon_max+1, lonid)
		! status = nf90_def_dim (ncid_to, "lat", 2*this.lat_max+1, latid)
		! status = nf90_def_dim (ncid_to, "time", all_time, time)
		! status = nf90_def_var (ncid_to, "water", NF90_DOUBLE, (/ lonid, latid, time/), Wid_to)
		! status = nf90_enddef (ncid_to)
		! if(status /= nf90_NoErr) print *, nf90_strerror(status)
		open(11,file=path3,access="direct",recl=(2*this.lat_max+1)*(2*this.lon_max+1))


		this.ncid = ncid;  this.ncid_gr = ncid_gr;  this.grid_id = grid_id(1); this.Wid = Wid(1)

	end subroutine



	subroutine scan_surf(this, time, surface_off)

		Class(printer) :: this
		integer(4), intent(in) :: time
		real(8), intent(out) :: surface_off(0:2*this.dim+1, 0:2*this.dim+1, 1:6)
		integer(4) x, y, face, ier, status, ncid, Wid, dim

		dim = this.dim;  ncid = this.ncid;  Wid = this.Wid

		status = nf90_get_var(ncid, Wid, surface_off(1:2*dim, 1:2*dim, 1:6),&
		 start = (/1, 1, 1, time/), count = (/2*dim, 2*dim, 6, 1/))
		if(status /= nf90_NoErr) print *, nf90_strerror(status), "scan_surf"
	end subroutine



	subroutine scan_grid(this, grid)

		Class(printer) :: this
		real(8), intent(out) :: grid(1:2, 0:2*this.dim+1, 0:2*this.dim+1, 1:6)
		integer(4) x, y, face, ier, dim, status, grid_id, ncid_gr

		dim = this.dim;  ncid_gr = this.ncid_gr;  grid_id = this.grid_id

		status = nf90_get_var(ncid_gr, grid_id, grid(1:2, 1:2*dim, 1:2*dim, 1:6),&
		 start = (/1, 1, 1, 1/), count = (/2, 2*dim, 2*dim, 6/))
		if(status /= nf90_NoErr) print *, nf90_strerror(status), "scan_grid"
	end subroutine



	subroutine print_surf(this, surface_to, time)
		Class(printer) :: this
		real(8), intent(in) :: surface_to(-this.lon_max:this.lon_max, -this.lat_max:this.lat_max)
		integer(4), intent(in) :: time
		integer(4) x, y, face, ier, status, Wid_to, ncid_to

		ncid_to = this.ncid_to;  Wid_to = this.Wid_to

		! status = nf90_put_var(ncid_to, Wid_to, surface_to(-this.lon_max:this.lon_max, -this.lat_max:this.lat_max),&
		!  start = (/1, 1, time/), count = (/2*this.lon_max+1, 2*this.lat_max+1, 1/))
		! if(status /= nf90_NoErr) print *, nf90_strerror(status)

		write(11, rec=time) real(surface_to(:,:),4)
	end subroutine



	subroutine deinit(this)
		Class(printer) :: this
		integer(4) :: ncid, ncid_to, ncid_gr, status
		ncid_gr = this.ncid_gr;  ncid = this.ncid; ncid_to = this.ncid_to

		status = nf90_close(ncid)
		status = nf90_close(ncid_gr)
		! status = nf90_close(ncid_to)
		close(11)
	end subroutine



end module