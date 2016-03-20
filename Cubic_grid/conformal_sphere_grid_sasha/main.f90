!!!!To know more look through Rancic Parser Mesinger article!!!! 
!!!!!!!!!Q.J.R. Meteorol. Soc. 1996, 122, pp. 959-982!!!!!!!!!!!

program main

	use grid_generator, Only: grid
	use data_analyzer, Only: data_an

	Type(grid) :: generator
	Type(data_an) :: analyzer

	call generator.conformal_cubed_sphere(50,50)
	call analyzer.generated()
	call analyzer.computated()

end program
