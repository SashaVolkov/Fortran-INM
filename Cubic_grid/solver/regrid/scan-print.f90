module scan_print

	use netcdf

	implicit none

	Private
	Public :: printer

	Type printer
		CONTAINS
		Procedure, Public :: init => init
		Procedure, Public :: to_print => to_print
		Procedure, Public :: deinit => deinit
	End Type

	CONTAINS


	subroutine init(this, dim, all_time, rescale, grid_type)

		Class(printer) :: this
		integer(4), intent(in) :: dim, all_time, rescale, grid_type

		integer(4) status, face
		character(40) istring
		character(80) path

		write(istring, *) dim


		if(grid_type == 0) then
			if(rescale == 0) then
				path = trim('datFiles/'//"surface_conf_simple_"//trim(adjustl(istring))//".nc")
			else if(rescale == 1) then
				path = trim('datFiles/'//"surface_conf_tan_"//trim(adjustl(istring))//".nc")
			else if(rescale == 2) then
				path = trim('datFiles/'//"surface_conf_exp_"//trim(adjustl(istring))//".nc")
			end if
		else if(grid_type == 1) then
				path = trim('datFiles/'//"surface_equiang_"//trim(adjustl(istring))//".nc")
		end if

		status = nf90_open (path = path,cmode = NF90_NOWRITE, ncid = ncid)


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

	end subroutine



	subroutine to_print(this, var, Wid, ncid, id)

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