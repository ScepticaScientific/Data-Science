module kronFunction

contains

	function kron(vector1, vector2)

	implicit none

	! Variables (partial)
	real vector1(:), vector2(:)
	real kron(size(vector1, 1), size(vector2, 1))
	integer i

	! Body
	kron(:, :) = 0.0

	!open(1, file = 'kron.txt')
	do i = 1, size(vector1, 1)
	   kron(i, :) = vector1(i) * vector2
	!   write (1, '(F20.15$)') (kron(i, j), j = 1, size(kron, 2) - 1)
	!   write (1, '(F20.15)') kron(i, size(kron, 2))
	end do
	!close(1)

	end function kron

end module kronFunction