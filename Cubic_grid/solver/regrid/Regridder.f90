program regrid

	use scan_print, Only: printer
	use grid_interp, Only: interp
	use sphere_geometry, Only: geometry
	use diagn, Only: diagnostic

	implicit none

	integer(4) dim, Tmax, speedup, rescale, grid_type, all_time, time
	real(8) dt, convert_time
	Type(printer) :: scan
	Type(interp) :: inter
	Type(geometry) :: geom
	Type(diagnostic) :: d


	open(9,file='../init')
		read(9, *) dim, Tmax, dt, speedup, rescale, grid_type
	close(9)

	all_time = Tmax/speedup + 1
	convert_time = dt

	call inter.init(dim)
	call scan.init(dim, all_time, convert_time, rescale, grid_type)
	call d.init(convert_time, grid_type, rescale)

	call scan.scan_grid(inter.latlon_c_off)
	call inter.weight_find(geom, scan)
	do time = 1, all_time
		call scan.scan_surf(time, inter.surface_off)
		call inter.interpolate()
		call scan.scan_precise(time, d.surface_precise)
		call d.L_norm((time-1)*speedup, inter.surface_to)
		call scan.print_surf(real(inter.surface_to,4), time)
	end do

	call inter.deinit()
	call scan.deinit()
	call d.deinit()

end program