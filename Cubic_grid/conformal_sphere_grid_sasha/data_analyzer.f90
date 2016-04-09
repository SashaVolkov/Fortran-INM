module data_analyzer

	use geometry
	implicit none

	Private
	Public :: data_an

	Type data_an
		CONTAINS
		Procedure :: generated => analyze_data_generation
		Procedure :: computated => analyze_data_computation
	End Type

	CONTAINS

	subroutine analyze_data_generation(this)

		Class(data_an) :: this

		integer x_points, y_points
		integer(4) i_face, i, j, status, k
		character*14 filename
		character istring
		real*8, dimension(:,:,:,:), allocatable :: r

		real*8 angle_cell(4), distance_cell(4), square


		open (20, file = "grid/parameters.dat")
		read(20,*) x_points, y_points
		close(20)
		allocate(r(1:3,1:6,-x_points:x_points, -y_points:y_points))

		do i_face = 1,6
			write(istring(1:1), '(i1.1)') i_face
			filename = "grid/face" // istring // ".dat" 
			open(20, file = filename)
			do i = -x_points, x_points
				do j = -y_points, y_points
					read(20,*) r(1,i_face,i,j), r(2,i_face,i,j), r(3,i_face,i,j)
				end do
			end do
			close(20)
		end do

		open(20, file = "analyze/parameters.dat")
		write(20,*) x_points, y_points
		close(20)

		open (20, file = "analyze/angle.dat")
		open (21, file = "analyze/distance.dat")
		open (22, file = "analyze/square.dat")
		do i = -x_points, x_points-1
			do j = -y_points, y_points-1
				status = cell_analyzer(r(:, 1, i, j), r(:, 1, i+1, j), r(:, 1, i, j+1), r(:, 1,i+1, j+1), &
				angle_cell(1:4), distance_cell(1:4), square)
				do k=1,4
					write(20,*) angle_cell(k)
					write(21,*) distance_cell(k)
				end do
				write(22,*) square
			end do
		end do

		close(20);close(21);close(22)

	end subroutine analyze_data_generation


	integer function cell_analyzer(r1,r2,r3,r4, angle_cell, distance_cell, square)
	! 
	! 1_____2
	! |     |
	! |_____|
	! 3     4

		real*8, dimension(1:3) :: r1,r2,r3,r4
		
		real*8 angle_cell(4), distance_cell(4), square

		angle_cell(1) = angle_r(r1,r2,r3)
		angle_cell(2) = angle_r(r2,r1,r4)
		angle_cell(3) = angle_r(r3,r1,r4)
		angle_cell(4) = angle_r(r4,r2,r3)

		distance_cell(1) = distance_r(r1,r2)
		distance_cell(2) = distance_r(r2,r4)
		distance_cell(3) = distance_r(r4,r3)
		distance_cell(4) = distance_r(r3,r1)

		square = 5d-1 * (distance_cell(1) *  distance_cell(4) * sin(angle_cell(1)) + &
										 distance_cell(2) *  distance_cell(3) * sin(angle_cell(4)))
	end function cell_analyzer


	subroutine analyze_data_computation(this)

	Class(data_an) :: this

	integer, parameter :: hist_points=500

	real*8, parameter :: pi = 314159265358979323d-17
	real*8, dimension(:), allocatable :: dat
	integer distribution(1:hist_points)
	real*8 max, min, max_val
	integer x_points, y_points
	integer i
	integer index

	open(20,file = "analyze/parameters.dat")
	read(20,*) x_points, y_points
	close(20)

	open(22, file="analyze/ratios.dat")
	call hist_generation(16*x_points*y_points, "analyze/angle.dat", "analyze/angle_distribution.dat", 180/pi)

	call hist_generation(16*x_points*y_points, "analyze/distance.dat", "analyze/distance_distribution.dat", 1d0)
	write(22,*) "distance ratio:", min/max

	call hist_generation(2*x_points*y_points, "analyze/square.dat", "analyze/square_distribution.dat", 1d0)
	write(22,*) "square ratio:", min/max
	close(22)

	CONTAINS
		subroutine hist_generation(N, filename_input, filename_output, x_coeff)
			real*8 x_coeff
			integer N
			character (len=*) filename_input, filename_output

			max=-1d20;min=1d20
			distribution = 0
			allocate(dat(1: N))
			open(20, file = filename_input)
			do i =1, N
				read(20,*) dat(i)
				if (dat(i)> max) max = dat(i)
				if (dat(i)< min) min = dat(i)
			end do
			close(20)

			do i = 1, N
				 index = int(hist_points * (dat(i)-min) /(max-min) + 5d-1)
				 distribution(index) = distribution(index) + 1
			end do

			max_val = maxval(distribution)

			open(21, file = filename_output)
			do i = 1, hist_points
				if (distribution(i)>0) then
					write(21,*) x_coeff * (min + (max-min) * i /dble(hist_points)), &
								distribution(i) !/ dble(max_val)
				end if
			end do
			close(21)
			deallocate(dat)
		end subroutine
	end subroutine


end module data_analyzer