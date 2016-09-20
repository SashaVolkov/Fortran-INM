program regrid

	use scan_print, Only: printer

	implicit none

	integer(4) dim, Tmax, speedup, rescale, grid_type, all_time
	Type(printer) :: scan


		open(9,file='../init')
			read(9, *) dim, Tmax, speedup, rescale, grid_type
		close(9)

		all_time = Tmax/speedup + 1

		call scan.init(dim, all_time, rescale, grid_type)



end program