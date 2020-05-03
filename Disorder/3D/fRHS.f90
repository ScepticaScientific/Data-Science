module fRHSFunction

contains

	! Third experiment for NMPDE (2010-2011) & First experiment for SIMULTECH 2011
	!function fRHS(omega1, omega2, phi1, phi2, xShift, time)

	! Nonlinear problem (for SIMULTECH 2011, the third experiment; for Int. J for Numerical Methods in Engineering, the fourth experiment; NOVA Publishers, the fourth experiment; WCE2013, first experiment)
	function fRHS(alpha, omega1, theta, kappa, the_am, the_add, myu_factor, time)

	! New nonlinear problem (for Int. J for Numerical Methods in Engineering and for Springer)
	!function fRHS(alpha, omega1, theta, kappa, phase, theta2, kappa2, phase2, gamma, the_am, the_add, myu_factor, time)

	implicit none

	! Variables (general)
	integer lr, lx, ly

	real, pointer :: x(:), y(:)

	real, pointer :: xs(:, :, :), ys(:, :, :)

    real, pointer :: mycosTabX(:), mysinTabX(:), mytanTabX(:)

	real, pointer :: cosphi(:, :), sinphi(:, :), tanphi(:, :)

	real, pointer :: myu(:, :, :)

	real R_Earth

	real pi 

	common /axes/ x, y
	common /axesmat/ xs, ys
	common /trigsX/ mycosTabX, mysinTabX, mytanTabX
	common /trigsXYmat/ cosphi, sinphi, tanphi
	common /sizes/ lr, lx, ly
	common /myu/ myu
	common /Radius/ R_Earth
	common /consts/ pi

	! Variables (partial)
	real fRHS(lr, lx, ly)
	real phi1, phi2, xShift
	real alpha
	real time
	real omega1, omega2
	real theta, theta2
	real kappa, kappa2
	real phase, phase2
	real gamma
	real the_am, the_add
	real myu_factor

	!real xi(lx, ly), dxidphi(lx, ly)					! Auxiliary variable for the new nonlinear problem

	integer j

	! Body

	! Third experiment for NMPDE (2010-2011) & First experiment for SIMULTECH 2011
	!fRHS = (-2.0 * myu * cos(time) ** 2.0 * (1.0 / spread(mycosTabX, 1, lx) * sin(omega2 * ys + phi2) ** 2.0 * (sin(omega1 * xs + phi1) ** 2.0 + 4.0 * omega1 * spread(x - xShift, 2, ly) * sin(omega1 * xs + phi1) * cos(omega1 * xs + phi1) + omega1 ** 2.0 * spread(x - xShift, 2, ly) ** 2.0 * cos(2.0 * (omega1 * xs + phi1))) + &
	!		omega2 * (spread(x - xShift, 2, ly) * sin(omega1 * xs + phi1)) ** 2.0 * (omega2 * spread(mycosTabX, 1, lx) * cos(2.0 * (omega2 * ys + phi2)) - 0.5 * spread(mysinTabX, 1, lx) * sin(2.0 * (omega2 * ys + phi2)))) / &
	!		R_Earth ** 2.0 / spread(mycosTabX, 1, lx)) - (spread(x - xShift, 2, ly) * sin(omega1 * xs + phi1)) ** 2.0 * sin(omega2 * ys + phi2) ** 2.0 * sin(2.0 * time)

	! Nonlinear problem (for SIMULTECH 2011, the third experiment; for Int. J for Numerical Methods in Engineering, the fourth experiment; NOVA Publishers, the fourth experiment; WCE2013, first experiment)
	fRHS = -the_am * myu_factor * theT(omega1, theta, kappa, the_am, the_add, time) ** (alpha - 1) * cos(time) ** 2.0 * (the_am * alpha * spread(cosphi, 3, lr) * cos(time) ** 2.0 * (omega1 ** 2.0 * cos(omega1 * xs + theta * cos(kappa * ys) * sin(time)) ** 2.0 + spread(sinphi, 3, lr) * sin(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * (kappa * theta * sin(time) * cos(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * sin(kappa * ys) * spread(cosphi, 3, lr) + sin(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * spread(sinphi, 3, lr))) + &
		   theT(omega1, theta, kappa, the_am, the_add, time) * sin(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * (spread(sinphi, 3, lr) ** 2.0 - spread(cosphi, 3, lr) ** 2.0 - omega1 ** 2.0) + & 
		   kappa * theta * sin(time) * spread(cosphi, 3, lr) * (alpha * the_am * cos(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * sin(kappa * ys) * spread(cosphi, 3, lr) * cos(time) ** 2.0 * (kappa * theta * sin(time) * cos(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * sin(kappa * ys) * spread(cosphi, 3, lr) + sin(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * spread(sinphi, 3, lr)) - theT(omega1, theta, kappa, the_am, the_add, time) * (kappa * theta * sin(time) * sin(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * sin(kappa * ys) ** 2.0 * spread(cosphi, 3, lr) + cos(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * (kappa * cos(kappa * ys) * spread(cosphi, 3, lr) - 3.0 * sin(kappa * ys) * spread(sinphi, 3, lr))))) / & 
		   R_Earth ** 2.0 / spread(spread(mycosTabX, 1, lx), 3, lr) + the_am * spread(cosphi, 3, lr) * (cos(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * cos(time) ** 3.0 * theta * cos(kappa * ys) - sin(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * sin(2.0 * time))
	!fRHS = 0.0

	! New nonlinear problem (for Int. J for Numerical Methods in Engineering and for Springer)
	!xi = omega1 * (xs - theta * tan(kappa * ys + phase)) + theta2 * tan(kappa2 * ys + phase2) * sin(gamma * time)
	!dxidphi = -omega1 * theta * kappa / cos(kappa * ys + phase) ** 2.0 + theta2 * kappa2 * sin(gamma * time) / cos(kappa2 * ys + phase2) ** 2.0

	!fRHS = the_am * cosphi * (cos(xi) * theta2 * tan(kappa2 * ys + phase2) * cos(gamma * time) * gamma * cos(time) ** 2.0 - sin(xi) * sin(2.0 * time)) - &
	!	   myu_factor * theT(omega1, theta, kappa, phase, theta2, kappa2, phase2, gamma, the_am, the_add, time) ** (alpha - 1.0) * ((alpha * (the_am * omega1 * cosphi * cos(time) ** 2.0 * cos(xi)) ** 2.0 - theT(omega1, theta, kappa, phase, theta2, kappa2, phase2, gamma, the_am, the_add, time) * the_am * omega1 ** 2.0 * cosphi * cos(time) ** 2.0 * sin(xi)) / R_Earth / spread(mycosTabX, 1, lx) + &
	!	   cosphi * (theT(omega1, theta, kappa, phase, theta2, kappa2, phase2, gamma, the_am, the_add, time) * (the_am * cos(time) ** 2.0 * (-2.0 * sinphi * cos(xi) * dxidphi - cosphi * sin(xi) * (1.0 + dxidphi ** 2.0) + 2.0 * cosphi * cos(xi) * (-omega1 * theta * kappa ** 2.0 * sin(kappa * ys + phase) / cos(kappa * ys + phase) ** 3.0 + theta2 * kappa2 ** 2.0 * sin(gamma * time) * sin(kappa2 * ys + phase2) / cos(kappa2 * ys + phase2) ** 3.0)) - tanphi * the_am * cos(time) ** 2.0 * (cosphi * cos(xi) * dxidphi - sinphi * sin(xi))) + alpha * (the_am * cos(time) ** 2.0 * (cosphi * cos(xi) * dxidphi - sinphi * sin(xi))) ** 2.0) / R_Earth) / R_Earth / spread(mycosTabX, 1, lx)

	end function fRHS

!!!!!!!

	! Third experiment for NMPDE (2010-2011) & First experiment for SIMULTECH 2011
	!function theT(omega1, omega2, phi1, phi2, xShift, time)

	! Nonlinear problem (for SIMULTECH 2011, the third experiment; for Int. J for Numerical Methods in Engineering, the fourth experiment; NOVA Publishers, the fourth experiment; WCE2013, first experiment)
	function theT(omega1, theta, kappa, the_am, the_add, time)

	! New nonlinear problem (for Int. J for Numerical Methods in Engineering and for Springer)
	!function theT(omega1, theta, kappa, phase, theta2, kappa2, phase2, gamma, the_am, the_add, time)

	implicit none

	! Variables (general)
	integer lr, lx, ly

	real, pointer :: x(:), y(:)

	real, pointer :: xs(:, :, :), ys(:, :, :)

	real, pointer :: cosphi(:, :), sinphi(:, :), tanphi(:, :)

	real pi 

	common /axes/ x, y
	common /axesmat/ xs, ys
	common /trigsXYmat/ cosphi, sinphi, tanphi
	common /sizes/ lr, lx, ly
	common /consts/ pi

	! Variables (partial)
	real theT(lr, lx, ly)
	real phi1, phi2, xShift
	real time
	real omega1, omega2
	real theta, theta2
	real kappa, kappa2
	real phase, phase2
	real gamma
	real the_am, the_add

	integer j

	! Third experiment for NMPDE (2010-2011)
	!theT = (spread(x - xShift, 2, ly) * sin(omega1 * xs + phi1)) ** 2.0 * sin(omega2 * ys + phi2) ** 2.0 * cos(time) ** 2.0 + 10.0

	! First experiment for SIMULTECH 2011
	!theT = (spread(x - xShift, 2, ly) * sin(omega1 * xs + phi1)) ** 2.0 * sin(omega2 * ys + phi2) ** 2.0 * cos(time) ** 2.0 + 1.0

	! Nonlinear problem (for SIMULTECH 2011, the third experiment; for Int. J for Numerical Methods in Engineering, the fourth experiment; NOVA Publishers, the fourth experiment; WCE2013, first experiment)
	theT = the_am * sin(omega1 * xs + theta * cos(kappa * ys) * sin(time)) * spread(cosphi, 3, lr) * cos(time) ** 2.0 + the_add
	!theT = 0.0

	! New nonlinear problem (for Int. J for Numerical Methods in Engineering and for Springer)
	!theT = the_am * sin(omega1 * (xs - theta * tan(kappa * ys + phase)) + theta2 * tan(kappa2 * ys + phase2) * sin(gamma * time)) * cosphi * cos(time) ** 2.0 + the_add

	end function theT


end module fRHSFunction