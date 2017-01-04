program regrid

	use scan_print, Only: printer
	use grid_interp, Only: interp
	use sphere_geometry, Only: geometry
	use diagn, Only: diagnostic

	implicit none

	integer(4) dim, Tmax, speedup, rescale, grid_type, all_time, time, step, lon_max, lat_max
	real(8) dt, convert_time, max
	Type(printer) :: scan
	Type(interp) :: inter
	Type(geometry) :: geom
	Type(diagnostic) :: d


	open(9,file='../init')
		read(9, *) dim, Tmax, dt, speedup, rescale, grid_type, step
	close(9)

	all_time = Tmax/speedup + 1
	convert_time = dt; lon_max = 360; lat_max = 180

	call inter.init(dim, lon_max, lat_max)
	call scan.init(dim, step, all_time, convert_time, rescale, grid_type, lon_max, lat_max)
	call d.init(dim, step, convert_time, grid_type, rescale, lon_max, lat_max)

	call scan.scan_grid(inter.latlon_c_off)
	call inter.weight_find(geom, scan)
	do time = 1, all_time
		call scan.scan_surf(time, inter.surface_off)
		call inter.interpolate(max)
		call scan.scan_precise(time, d.surface_precise)
		call d.L_norm((time-1)*speedup, inter.surface_to, max)
		call scan.print_surf(inter.surface_to, d.surface_precise, time)
	end do

	call inter.deinit()
	call scan.deinit()
	call d.deinit()

end program