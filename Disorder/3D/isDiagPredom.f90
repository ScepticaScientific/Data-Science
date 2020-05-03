module isDiagPredomFunction

contains

	function isDiagPredom(matrix)

	implicit none

	! Variables (partial)
	real matrix(:, :)
	logical isDiagPredom
	integer i

	! Body
	do i = 2, size(matrix, 1) - 1
	   if (sum(abs([matrix(i, 1 : i - 1), matrix(i, i + 1 : size(matrix, 2))])) >= abs(matrix(i, i))) then
	      isDiagPredom = .false.
		  return
	   end if
    end do

	isDiagPredom = .true.

	end function isDiagPredom

end module isDiagPredomFunction