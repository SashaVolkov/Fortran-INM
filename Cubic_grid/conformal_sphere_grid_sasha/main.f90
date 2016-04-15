!!!!To know more look through Rancic Parser Mesinger article!!!! 
!!!!!!!!!Q.J.R. Meteorol. Soc. 1996, 122, pp. 959-982!!!!!!!!!!!

program main

	use grid_generator, Only: grid
	use data_analyzer, Only: data_an

	real(8) t(4)
	integer(4) dimention

	Type(grid) :: generator
	Type(data_an) :: analyzer

	dimention = 100

	call cpu_time(t(1))
	call generator.conformal_cubed_sphere(dimention,dimention)
	call cpu_time(t(2))
	print '("Time of generation = ", f6.3, " sec")', t(2) - t(1)

	call analyzer.generated()
	call cpu_time(t(3))
	print '("Time of 1 analyze step = ", f6.3, " sec")', t(2) - t(1)

	call analyzer.computated()
	call cpu_time(t(4))
	print '("Time of 2 analyze step = ", f6.3, " sec")', t(2) - t(1)

end program
