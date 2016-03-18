module simple_rotations
	real(8), dimension(1:3,1:3) :: Rotation_xy = reshape(source = &
					 (/0, 1, 0, &
						-1, 0, 0, &
						 0, 0, 1 /), shape=(/3,3/))
	
	real(8), dimension(1:3,1:3) :: Rotation_yz = reshape(source = &
			 (/ (/ 1, 0, 0 /), &
				 (/ 0, 0, 1 /), &
				 (/ 0,-1, 0 /)  /), shape = (/3,3/))

	real(8), dimension(1:3,1:3) ::  Rotation_mir = reshape(source = &
			 (/ (/ 1, 0, 0 /), &
			 (/ 0,-1, 0 /), &
			 (/ 0, 0, 1 /)  /), shape = (/3,3/))
end module simple_rotations