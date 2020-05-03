module GetMassSubroutine

contains

	subroutine GetMass(mass, energy)

	implicit none

	! Variables (general)
	!real, pointer :: a(:, :, :), rhs(:, :, :)
	real, pointer :: a(:, :, :)

	real, pointer :: mycosTabX(:), mysinTabX(:), mytanTabX(:)

	real dr, dx, dy
	integer lr, lx, ly

	real tau

	real R_Earth

	!common /arhs/ a, rhs
	common /a/ a
	common /trigsX/ mycosTabX, mysinTabX, mytanTabX
	common /steps/ dr, dx, dy
	common /sizes/ lr, lx, ly
	common /tau/ tau
	common /Radius/ R_Earth

	! Variables (partial)
	real mass, energy

	integer ij_r, j

	! Body
	mass = 0.0
	energy = 0.0

	do ij_r = 1, lr
	   do j = 1, ly
          mass = mass + mycosTabX(j) * sum(a(1 : lx - 2, j, ij_r))
	      energy = energy + mycosTabX(j) * sum(a(1 : lx - 2, j, ij_r) ** 2.0)
	   end do
	end do

	mass = mass * R_Earth * R_Earth * dr * dx * dy
	energy = sqrt(energy * R_Earth * R_Earth * dr * dx * dy)

	end subroutine GetMass

end module GetMassSubroutine