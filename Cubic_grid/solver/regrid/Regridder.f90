program regrid

	use scan_print, Only: printer
	use grid_interp, Only: interp
	use sphere_geometry, Only: geometry

	implicit none

	integer(4) dim, Tmax, speedup, rescale, grid_type, all_time, time
	Type(printer) :: scan
	Type(interp) :: inter
	Type(geometry) :: geom


	open(9,file='../init')
		read(9, *) dim, Tmax, speedup, rescale, grid_type
	close(9)

	all_time = Tmax/speedup + 1

	call inter.init(dim)
	call scan.init(dim, all_time, rescale, grid_type)
	call scan.scan_grid(inter.latlon_c_off)
	call inter.weight_find(geom)
	do time = 1, all_time
		call scan.scan_surf(time, inter.surface_off)
		call inter.interpolate()
		call scan.print_surf(inter.surface_to, time)
	end do

	call inter.deinit()
	call scan.deinit()

end program