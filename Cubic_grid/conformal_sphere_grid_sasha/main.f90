program main
use grid_generator
use data_analyzer

status = conformal_cubed_sphere_grid_generation(50,50)
!status = cubed_sphere_grid_generation(50,50)
status = analyze_data_generation()
ostatus = analyze_data_computation()

end program
