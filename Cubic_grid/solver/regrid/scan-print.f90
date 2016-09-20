module scan_print

	use netcdf

	implicit none

	Private
	Public :: printer

	integer(4) :: dim

	Type printer
		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: scan_surf => scan_surf
		Procedure, Public :: scan_grid => scan_grid
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS



	subroutine init(this, dim, all_time, rescale, grid_type)

		Class(printer) :: this
		integer(4), intent(in) :: dim, all_time, rescale, grid_type

		integer(4) status, face
		character(40) istring
		character(80) path1, path2

		this.dim = dim

		write(istring, *) 2*dim

		if(grid_type == 0) then
			if(rescale == 0) then
				path1 = trim('../datFiles/'//"surface_conf_simple_C"//trim(adjustl(istring))//".nc")
				path2 = trim('../datFiles/'//"grid_conf_simple_C"//trim(adjustl(istring))//".nc")
			else if(rescale == 1) then
				path1 = trim('../datFiles/'//"surface_conf_tan_C"//trim(adjustl(istring))//".nc")
				path2 = trim('../datFiles/'//"grid_conf_tan_C"//trim(adjustl(istring))//".nc")
			else if(rescale == 2) then
				path1 = trim('../datFiles/'//"surface_conf_exp_C"//trim(adjustl(istring))//".nc")
				path2 = trim('../datFiles/'//"grid_conf_exp_C"//trim(adjustl(istring))//".nc")
			end if
		else if(grid_type == 1) then
				path1 = trim('../datFiles/'//"surface_equiang_C"//trim(adjustl(istring))//".nc")
				path2 = trim('../datFiles/'//"grid_equiang_C"//trim(adjustl(istring))//".nc")
		end if


		status = nf90_open (path = path1,cmode = NF90_NOWRITE, ncid = ncid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		status = nf90_def_dim (ncid, "x", 2*dim, xid)
		status = nf90_def_dim (ncid, "y", 2*dim, yid)
		status = nf90_def_dim (ncid, "face", 6, faceid)
		status = nf90_def_dim (ncid, "time", all_time, time)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)

		status = nf90_def_var (ncid, "water", NF90_DOUBLE, (/ xid, yid, faceid, time/), Wid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)
		status = nf90_enddef (ncid)
		if(status /= nf90_NoErr) print *, nf90_strerror(status)


		status = nf90_open (path = path2,cmode = NF90_NOWRITE, ncid = ncid_gr)
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


	end subroutine



	subroutine scan_surf(this, Wid, ncid)

		Class(printer) :: this
		integer(4), intent(in) :: time, speedup, Wid, ncid
		real(8), intent(in) :: h_height(1:2*dim, 1:2*dim, 1:6)
		integer(4) x, y, face, ier
		integer(4) status, t, ns_y, ns_x, nf_y, nf_x, Ysize, Xsize

		t = 1+time/speedup

		status = nf90_get_var(ncid, Wid, h_height(1:2*dim, 1:2*dim, 1:6),&
		 start = (/1, 1, 1/), count = (/2*dim, 2*dim, 6/))

		if(status /= nf90_NoErr) print *, nf90_strerror(status)
	end subroutine



	subroutine scan_grid(this, grid, grid_id, ncid_gr)

		Class(printer) :: this
		Class(g_var) :: grid
		integer(4), intent(in) :: grid_id, ncid_gr
		integer(4) x, y, face, ier, dim, status

		dim = this.dim

		status = nf90_get_var(ncid_gr, grid_id, grid.latlon_c(1:2, 1:2*dim, 1:2*dim, 1:6),&
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