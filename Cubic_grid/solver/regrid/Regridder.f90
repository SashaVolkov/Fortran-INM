program regrid

	use scan_print, Only: printer
	use grid_interp, Only: interp
	use sphere_geometry, Only: geometry
	use diagn, Only: diagnostic
	use interp_prec_to_cube, Only: prec_to_cube

	implicit none

	integer(4) dim, Tmax, speedup, rescale, grid_type, all_time, time, step, lon_max, lat_max, nc_or_dat
	real(8) dt, convert_time, max
	Type(printer) :: scan
	Type(interp) :: inter
	Type(geometry) :: geom
	Type(diagnostic) :: d
	Type(prec_to_cube) :: to_cube


	open(9,file='../init')
		read(9, *) dim, Tmax, dt, speedup, rescale, grid_type, step
	close(9)

	all_time = Tmax/speedup + 1
	convert_time = dt; lon_max = 360; lat_max = 180
	nc_or_dat = 1

	call inter.init(dim, step, lon_max, lat_max)
	call to_cube.init(dim, step, lon_max, lat_max)
	call scan.init(dim, step, all_time, convert_time, rescale, grid_type, lon_max, lat_max, nc_or_dat)
	call d.init(dim, step, convert_time, grid_type, rescale, lon_max, lat_max)

	call scan.scan_grid(inter.latlon_cubic)
	call inter.weight_find(geom, scan)
	to_cube.latlon_cubic = inter.latlon_cubic
	call to_cube.weight_find()
	do time = 1, all_time
		call scan.scan_surf(time, inter.surface_off)
		call scan.scan_precise(time, d.surface_precise)
		call inter.interpolate(max)
		call to_cube.interpolate(d.surface_precise)
! 		call d.L_norm_ll((time-1)*speedup, inter.surface_to, max)
		call d.L_norm_c((time-1)*speedup, inter.surface_off, to_cube.precise_cube, max)
		call scan.print_surf(inter.surface_to, d.surface_precise, time)
		call scan.print_surf_prec_cube(time, to_cube.precise_cube, inter.surface_off)
	end do

	call to_cube.deinit()
	call inter.deinit()
	call scan.deinit()
	call d.deinit()

end program