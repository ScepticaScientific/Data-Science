program Diffusion3DSphere

! To do: Use 'MYU_FACTOR_IS_ARRAY' if myu is a function of (lambda, phi) or use 'NOTMYU_FACTOR_IS_ARRAY' if myu is a constant
!DEC$ DEFINE NOTMYU_FACTOR_IS_ARRAY

! To do: Use 'SOLUTION_ENERGY' if you need to compute the solution's L2-norm or use 'ERROR_ENERGY' if you need to compute the solution relative error's L2-norm (i.e. ||(numerics - analytics) / analytics||_L2)
!DEC$ DEFINE SOLUTION_ENERGY

!!! IMPORTANT NOTE ABOUT NONLINEAR DIFFUSION (alpha .ne. 0) !!!
! THE PARAMETER a AT THE INITIAL TIME MOMEMT t = 0 SHOULD NOT EXCEED 1.0.
! IF IT DOES EXCEED, THEN THE NONLINEAR TERM myu_factor * a ** alpha MAY DRASTICALLY INCREASE THE SOLUTION (ESPECIALLY WHEN MODELLING COMBUSTION, 
! i.e. f ~ a ** beta AT beta > 1), AND THE FINITE DIFFERENCE SCHEMES, IN ORDER TO STAY ACCURATE, WILL REQUIRE THE TIME STEP tau TO BE VERY SMALL ---
! OTHERWISE THE GRID NUMBERS cflR, cflX AND cflY MAY BE ESSENTIALLY GREATER THAN 1.0.

use linear_operators

use kronFunction
use GetMassSubroutine
use fRHSFunction

implicit none

!DEC$ IF DEFINED (MPI)
include 'mpif.h'
!!! IMPORTANT NOTE: IF THE MPI IS USED THEN THE NUMBER OF PROCESSES MUST BE SUCH THAT THEY CONTAIN THE SAME NUMBER OF GRID LINES PER PROCESS
!!! IN ALL THREE SUBROUTINES 'DoInR()', 'DoInX()' and 'DoInY()' (SEE THE COMMENTS THERE TOO). OTHERWISE THE SUBROUTINE 'mpi_gather()' WILL 
!!! PROVIDE INCORRECT GATHERING RESULTS, SEE THE DOCUMENTATION ON 'mpi_gather()'
!DEC$ ENDIF

! Variables (general)
!real, pointer :: a(:, :, :), rhs(:, :, :)
real, pointer :: a(:, :, :)
!!!real, pointer :: rhs(:, :, :)

real, pointer :: xy_r(:), x(:), y(:)

!!!real, pointer :: xs(:, :, :), ys(:, :, :)

real, pointer :: my_r2(:)
real, pointer :: mycosTabX(:), mysinTabX(:), mytanTabX(:)
real, pointer :: mycosTabY(:), mysinTabY(:), mytanTabY(:)

!!!real, pointer :: cosphi(:, :), sinphi(:, :), tanphi(:, :)

real dr, dx, dy
integer lr, lx, ly
real llr, llx, lly					! These are Lr, Lx and Ly, respectively, in MatLab

real, pointer :: yNew(:)
integer lyNew
real llyNew							! This is LyNew in MatLab

real, pointer :: myu(:, :, :)

real R_Earth

real pi

real mass_dummy, energy_aux
real alpha

real, pointer :: Ar(:, :), bra(:), ipivr(:)
real, pointer :: Ax(:, :), bxa(:), ipivx(:)
real, pointer :: Ay(:, :), bya(:), ipivy(:)

!DEC$ IF DEFINED (MYU_FACTOR_IS_ARRAY)
real, pointer :: myu_factor(:, :, :)
!DEC$ ENDIF
!DEC$ IF .NOT. DEFINED (MYU_FACTOR_IS_ARRAY)
real myu_factor
!DEC$ ENDIF

integer, pointer :: indexMainDiagR(:), indexSuperDiagR(:), indexSubDiagR(:)
!integer, pointer :: indexMainDiagX(:), indexSuperDiagX(:), indexSuperDiagX2(:), indexSuperDiagX3(:), indexSuperDiagX4(:), indexSubDiagX(:), indexSubDiagX2(:), indexSubDiagX3(:), indexSubDiagX4(:)
!integer, pointer :: indexMainDiagY(:), indexSuperDiagY(:), indexSuperDiagY2(:), indexSuperDiagY3(:), indexSuperDiagY4(:), indexSubDiagY(:), indexSubDiagY2(:), indexSubDiagY3(:), indexSubDiagY4(:)
integer, pointer :: indexMainDiagX(:), indexSuperDiagX(:), indexSubDiagX(:)
integer, pointer :: indexMainDiagY(:), indexSuperDiagY(:), indexSubDiagY(:)

real, pointer :: workingArrayR(:), workingArrayX(:), workingArrayY(:)

real cflR, cflX, cflY

!DEC$ IF DEFINED (MPI)
integer ierr, processID, NumberOfProcesses
!DEC$ ENDIF

!common /arhs/ a, rhs
common /a/ a
!!!common /rhs/ rhs
common /axes/ xy_r, x, y
!!!common /axesmat/ xs, ys
common /nontrigR/ my_r2
common /trigsX/ mycosTabX, mysinTabX, mytanTabX
common /trigsY/ mycosTabY, mysinTabY, mytanTabY
!!!common /trigsXYmat/ cosphi, sinphi, tanphi
common /steps/ dr, dx, dy
common /limits/ llr, llx, lly
common /sizes/ lr, lx, ly
common /tau/ tau
common /myu/ myu
common /Radius/ R_Earth
!common /specials/ ApprOrder
common /consts/ pi
common /matvecR/ Ar, bra, ipivr
common /diagsR/ indexMainDiagR, indexSuperDiagR, indexSubDiagR, workingArrayR
common /matvecX/ Ax, bxa, ipivx
!common /diagsX/ indexMainDiagX, indexSuperDiagX, indexSuperDiagX2, indexSuperDiagX3, indexSuperDiagX4, indexSubDiagX, indexSubDiagX2, indexSubDiagX3, indexSubDiagX4, workingArrayX
common /diagsX/ indexMainDiagX, indexSuperDiagX, indexSubDiagX, workingArrayX
common /matvecY/ Ay, bya, ipivy
!common /diagsY/ indexMainDiagY, indexSuperDiagY, indexSuperDiagY2, indexSuperDiagY3, indexSuperDiagY4, indexSubDiagY, indexSubDiagY2, indexSubDiagY3, indexSubDiagY4, workingArrayY
common /diagsY/ indexMainDiagY, indexSuperDiagY, indexSubDiagY, workingArrayY
common /CFL/ cflR, cflX, cflY

!DEC$ IF DEFINED (MPI)
common /mpiData/ ierr, processID, NumberOfProcesses
!DEC$ ENDIF

! Variables (partial)
real T, tau
real a_ampl

integer Nx, Ny, Nr

integer fieldNo
!logical IsDepth

integer currTime
real, pointer :: mass(:), energy(:)

real, pointer :: analytic_solution(:, :, :)

integer outputStep, outputToDiskStep
integer i, j, ij_r
integer res
logical IsSources

integer cx

! Body
!DEC$ IF DEFINED (MPI)
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD, processID, ierr)
call mpi_comm_size(MPI_COMM_WORLD, NumberOfProcesses, ierr)

write (*, *) 'NoP: ', NumberOfProcesses, '; PID: ', processID, '; ierr: ', ierr

!DEC$ ENDIF

pi = 3.141592653589793238D+000

! Definitions
R_Earth = 1.0 / (2.0 * pi)     ! Earth's radius

!write (*, *) 'Test.f90(): PID:', processID, '; T: ', T

!DEC$ IF DEFINED (MPI)
if (processID .eq. 0) then
!DEC$ ENDIF

open (1, file = 'init.txt')
read (1, *) T, tau, outputStep, outputToDiskStep
read (1, *) Nx, Ny, Nr
read (1, *) IsSources
close (1)

!DEC$ IF DEFINED (MPI)
end if

!write (*, *) 'before mpi_bcast(): PID:', processID, '; T: ', T, '; Nr: ', Nr

call mpi_bcast(T, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

!write (*, *) 'Test.f90(): PID:', processID, '; T: ', T, '; Nr: ', Nr, '; ierr: ', ierr

call mpi_bcast(tau, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

call mpi_bcast(Nx, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
call mpi_bcast(Ny, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
call mpi_bcast(Nr, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)

call mpi_bcast(IsSources, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

!DEC$ ENDIF

!write (*, *) 'after mpi_bcast(): PID:', processID, '; T: ', T, '; Nr: ', Nr, '; ierr: ', ierr

!llr = 1.0 / 3600.0; 		dr = llr / (Nr - 1)				! '1.0 / 3600.0' is approx. 10 km
llr = R_Earth; 				dr = llr / (Nr - 1)
llx = 2.0 * pi;  			dx = llx / (Nx - 1)
lly = pi / 2.0;  			dy = 2.0 * lly / (Ny - 1)

!outputStep = max(int(5 * 5.0e-4 / tau), 50)

allocate(xy_r(Nr));											lr = size(xy_r, 1)
xy_r = R_Earth + [0 : lr - 1] * dr							! Earth's radius is '1.0 / (2.0 * pi)'

allocate(x(int((llx + dx) / dx + 1)));						lx = size(x, 1)
x = [0 : lx - 1] * dx + dx / 2.0

allocate(y(int((2.0 * lly - dy) / dy + 1)));				ly = size(y, 1)
y = -lly + dy / 2.0 + [0 : ly - 1] * dy

!
!!!allocate(xs(lr, lx, ly), ys(lr, lx, ly))
!!!xs = spread(spread(x, 1, ly), 3, lr)
!!!ys = reshape(spread(reshape(spread(y, 1, lx), (/1, lx * ly/)), 1, lr), (/lr, lx, ly/))

!!!allocate(cosphi(lx, ly), sinphi(lx, ly), tanphi(lx, ly))
!!!cosphi = cos(ys(1, :, :))
!!!sinphi = sin(ys(1, :, :))
!!!tanphi = tan(ys(1, :, :))

allocate(my_r2(lr))
my_r2 = xy_r ** 2.0

allocate(mycosTabX(ly), mysinTabX(ly), mytanTabX(ly))
mycosTabX = cos(y)
mysinTabX = sin(y)
mytanTabX = tan(y)

! The parameters "LyNew" and "yNew" must be the same as those used in the function "SwapXY()"
llyNew = 2.0 * pi
allocate(yNew(int((llyNew + dy) / dy + 1)));				lyNew = size(yNew, 1)
yNew = -lly + dy / 2.0 + [0 : lyNew - 1] * dy

allocate(mycosTabY(lyNew), mysinTabY(lyNew), mytanTabY(lyNew))
mycosTabY = abs(cos(yNew))
mysinTabY = sin(yNew)
mytanTabY = tan(yNew)

deallocate(yNew)

!
allocate(myu(lx, ly, lr))

!DEC$ IF DEFINED (MYU_FACTOR_IS_ARRAY)
allocate(myu_factor(lx, ly, lr))
!DEC$ ENDIF

!
allocate(a(lx, ly, lr))
!!!allocate(rhs(lx, ly, lr), analytic_solution(lx, ly, lr))

!!!do ij_r = 1, size(rhs, 3)
!!!   rhs(:, :, ij_r) = 0.0 * kron(sin(x) ** 2.0, sin(y) ** 2.0)
!!!end do

!DEC$ IF DEFINED (MPI)
if (processID .eq. 0) then
!DEC$ ENDIF

! Saving the axes
open (1, file = 'x.tx2')
write (1, '(D20.10)') x
close (1)

open (1, file = 'y.tx2')
write (1, '(D20.10)') y
close (1)

!DEC$ IF DEFINED (MPI)
end if
!DEC$ ENDIF

deallocate(xy_r, x, y)

!DEC$ IF DEFINED (MPI)
if (processID .eq. 0) then
!DEC$ ENDIF

! Sources.
!!!if (IsSources .eq. 1) then
!!!   open (1, file = 'SourcesInit.txt')
!!!   do ij_r = 1, size(rhs, 3)
!!!      do i = 1, size(rhs, 1)
!!!	     read (1, *) (rhs(i, j, ij_r), j = 1, size(rhs, 2))
!!!      end do 
!!!   end do
!!!   close (1)
!!!end if

!DEC$ IF DEFINED (MPI)
end if

!!!call mpi_bcast(rhs, product(shape(rhs)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
!DEC$ ENDIF

!DEC$ IF DEFINED (MPI)
if (processID .eq. 0) then
!DEC$ ENDIF

open (1, file = 'aInit.txt')
do ij_r = 1, size(a, 3) 
   do i = 1, size(a, 1)
      read (1, *) (a(i, j, ij_r), j = 1, size(a, 2))
   end do
end do 
close (1)

open (1, file = 'myuInit.txt')
do ij_r = 1, size(myu, 3) 
   do i = 1, size(myu, 1)
      read (1, *) (myu(i, j, ij_r), j = 1, size(myu, 2))
   end do
end do
close (1)

!DEC$ IF DEFINED (MYU_FACTOR_IS_ARRAY)
myu_factor = myu
!DEC$ ENDIF

!DEC$ IF DEFINED (MPI)
end if

call mpi_bcast(a, product(shape(a)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
call mpi_bcast(myu, product(shape(myu)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

!DEC$ IF DEFINED (MYU_FACTOR_IS_ARRAY)
call mpi_bcast(myu_factor, product(shape(myu_factor)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
!DEC$ ENDIF

!write (*, *) processID, un(lx, ly), product(shape(un)), shape(un), vn(lx, ly), a(lr, lx, ly), U(lx, ly), V(lx, ly)
!DEC$ ENDIF

! These allocations and assignements are equal to those performed in MatLab
allocate(Ar(lr, lr), bra(lr), ipivr(lr))

allocate(indexSuperDiagR(size(Ar, 2) - 1 - 1))
indexSuperDiagR = (/(size(Ar, 1) * (j - 1) + j - 1, j = 3, size(Ar, 2))/)

allocate(indexMainDiagR(size(Ar, 2) - 1 - 1))
indexMainDiagR = (/(size(Ar, 1) * (j - 1) + j, j = 2, size(Ar, 2) - 1)/)

allocate(indexSubDiagR(size(Ar, 2) - 1 - 1))
indexSubDiagR = (/(size(Ar, 1) * (j - 1) + j + 1, j = 1, size(Ar, 2) - 2)/)

allocate(Ax(lx, lx), bxa(lx), ipivx(lx))

!allocate(indexSuperDiagX4(size(Ax, 2) - 1 - 1))
!indexSuperDiagX4 = (/(size(Ax, 1) * (j - 1) + j - 4, j = 6, size(Ax, 2)), size(Ax, 1) * (4 - 1) - 3, size(Ax, 1) * (5 - 1) - 2, size(Ax, 1) * (6 - 1) - 1/)

!allocate(indexSuperDiagX3(size(Ax, 2) - 1 - 1))
!indexSuperDiagX3 = (/(size(Ax, 1) * (j - 1) + j - 3, j = 5, size(Ax, 2)), size(Ax, 1) * (4 - 1) - 2, size(Ax, 1) * (5 - 1) - 1/)

!allocate(indexSuperDiagX2(size(Ax, 2) - 1 - 1))
!indexSuperDiagX2 = (/(size(Ax, 1) * (j - 1) + j - 2, j = 4, size(Ax, 2)), size(Ax, 1) * (4 - 1) - 1/)

allocate(indexSuperDiagX(size(Ax, 2) - 1 - 1))
indexSuperDiagX = (/(size(Ax, 1) * (j - 1) + j - 1, j = 3, size(Ax, 2))/)

allocate(indexMainDiagX(size(Ax, 2) - 1 - 1))
indexMainDiagX = (/(size(Ax, 1) * (j - 1) + j, j = 2, size(Ax, 2) - 1)/)

allocate(indexSubDiagX(size(Ax, 2) - 1 - 1))
indexSubDiagX = (/(size(Ax, 1) * (j - 1) + j + 1, j = 1, size(Ax, 2) - 2)/)

!allocate(indexSubDiagX2(size(Ax, 2) - 1 - 1))
!indexSubDiagX2 = (/size(Ax, 1) * (size(Ax, 2) - 2 - 1) + 2, (size(Ax, 1) * (j - 1) + j + 2, j = 1, size(Ax, 2) - 3)/)

!allocate(indexSubDiagX3(size(Ax, 2) - 1 - 1))
!indexSubDiagX3 = (/size(Ax, 1) * (size(Ax, 2) - 2 - 2) + 2, size(Ax, 1) * (size(Ax, 2) - 2 - 1) + 3, (size(Ax, 1) * (j - 1) + j + 3, j = 1, size(Ax, 2) - 4)/)

!allocate(indexSubDiagX4(size(Ax, 2) - 1 - 1))
!indexSubDiagX4 = (/size(Ax, 1) * (size(Ax, 2) - 2 - 3) + 2, size(Ax, 1) * (size(Ax, 2) - 2 - 2) + 3, size(Ax, 1) * (size(Ax, 2) - 2 - 1) + 4, (size(Ax, 1) * (j - 1) + j + 4, j = 1, size(Ax, 2) - 5)/)

allocate(Ay(lyNew, lyNew), bya(lyNew), ipivy(lyNew))

!allocate(indexSuperDiagY4(size(Ay, 2) - 1 - 1))
!indexSuperDiagY4 = (/(size(Ay, 1) * (j - 1) + j - 4, j = 6, size(Ay, 2)), size(Ay, 1) * (4 - 1) - 3, size(Ay, 1) * (5 - 1) - 2, size(Ay, 1) * (6 - 1) - 1/)

!allocate(indexSuperDiagY3(size(Ay, 2) - 1 - 1))
!indexSuperDiagY3 = (/(size(Ay, 1) * (j - 1) + j - 3, j = 5, size(Ay, 2)), size(Ay, 1) * (4 - 1) - 2, size(Ay, 1) * (5 - 1) - 1/)

!allocate(indexSuperDiagY2(size(Ay, 2) - 1 - 1))
!indexSuperDiagY2 = (/(size(Ay, 1) * (j - 1) + j - 2, j = 4, size(Ay, 2)), size(Ay, 1) * (4 - 1) - 1/)

allocate(indexSuperDiagY(size(Ay, 2) - 1 - 1))
indexSuperDiagY = (/(size(Ay, 1) * (j - 1) + j - 1, j = 3, size(Ay, 2))/)

allocate(indexMainDiagY(size(Ay, 2) - 1 - 1))
indexMainDiagY = (/(size(Ay, 1) * (j - 1) + j, j = 2, size(Ay, 2) - 1)/)

allocate(indexSubDiagY(size(Ay, 2) - 1 - 1))
indexSubDiagY = (/(size(Ay, 1) * (j - 1) + j + 1, j = 1, size(Ay, 2) - 2)/)

!allocate(indexSubDiagY2(size(Ay, 2) - 1 - 1))
!indexSubDiagY2 = (/size(Ay, 1) * (size(Ay, 2) - 2 - 1) + 2, (size(Ay, 1) * (j - 1) + j + 2, j = 1, size(Ay, 2) - 3)/)

!allocate(indexSubDiagY3(size(Ay, 2) - 1 - 1))
!indexSubDiagY3 = (/size(Ay, 1) * (size(Ay, 2) - 2 - 2) + 2, size(Ay, 1) * (size(Ay, 2) - 2 - 1) + 3, (size(Ay, 1) * (j - 1) + j + 3, j = 1, size(Ay, 2) - 4)/)

!allocate(indexSubDiagY4(size(Ay, 2) - 1 - 1))
!indexSubDiagY4 = (/size(Ay, 1) * (size(Ay, 2) - 2 - 3) + 2, size(Ay, 1) * (size(Ay, 2) - 2 - 2) + 3, size(Ay, 1) * (size(Ay, 2) - 2 - 1) + 4, (size(Ay, 1) * (j - 1) + j + 4, j = 1, size(Ay, 2) - 5)/)

allocate(workingArrayR(product(shape(Ar))), workingArrayX(product(shape(Ax))), workingArrayY(product(shape(Ay)))) 

allocate(mass(int(T / tau)), energy(int(T / tau)))

res = 1

cx = 0

!DEC$ IF DEFINED (MPI)
if (processID .eq. 0) then
!DEC$ ENDIF

call GetMass(mass(1), energy(1))
call SaveToDisk(cx)

!DEC$ IF DEFINED (MPI)
end if

call mpi_barrier(MPI_COMM_WORLD, ierr)
!DEC$ ENDIF

alpha = 1.0
!DEC$ IF .NOT. DEFINED (MYU_FACTOR_IS_ARRAY)
myu_factor = 1.0e-03
!DEC$ ENDIF

!write (*, *) 'Start'

! Computing
do currTime = 1, int(T / tau)
   cflR = -1.0D+200
   cflX = -1.0D+200
   cflY = -1.0D+200

   myu = myu_factor * a ** alpha

   Ar = 0.0
   bra = 0.0

   ! Computing in "r"
   call DoInR(tau / 2.0, res)

   !write (*, *) 'Press any key after R...'
   !pause

   Ax = 0.0
   bxa = 0.0

   ! Computing in "x"
   call DoInX(tau / 2.0, res)

   !write (*, *) 'Press any key...'
   !pause

   if (res .eq. 0) then
      pause
      stop
   end if

   call SwapXY()

   Ay = 0.0
   bya = 0.0

   ! Computing in "y"
   call DoInY(tau / 2.0, res)

   if (res .eq. 0) then
      pause
      stop
   end if
    
   call SwapYX()

   ! Sources
   if (IsSources .eq. 1) then
!!!	  a = a + tau * rhs
   else if (IsSources .eq. 2) then
!DEC$ IF DEFINED (ERROR_ENERGY)
      ! Third experiment for NMPDE (2010-2011)
	  !a = a + tau * fRHS(1.0, 2.0, 0.0, -pi / 3.0, pi, currTime * tau - tau / 2.0)

      ! First experiment for SIMULTECH 2011
	  !a = a + tau * fRHS(1.0, 1.0, 0.0, 0.0, pi, currTime * tau - tau / 2.0)

	  ! Nonlinear problem (for SIMULTECH 2011, the third experiment; for Int. J for Numerical Methods in Engineering, the fourth experiment; NOVA Publishers, the fourth experiment; WCE2013, first experiment)
	  a = a + tau * fRHS(alpha, 9.0, 5.0 * alpha, 3.0, -2.5, 50.0, myu_factor, currTime * tau - tau / 2.0)

	  ! New nonlinear problem (for Int. J for Numerical Methods in Engineering and for Springer)
	  !a = a + tau * fRHS(alpha, 7.0, 1.0, 0.5, 0.0 * pi / 4.0, 15.0, 0.7, 0.0 * pi / 4.0, 4.0, 10.0, 100.0, myu_factor, currTime * tau - tau / 2.0)
!DEC$ ENDIF
      ! NOVA Publishers, third experiment (combustion); Applied Mathematics and Computation, fifth experiment (combustion); WCE2013, temperature wave; WCE2013, combustion
	  !a = a + tau * 1.0e+01 * (a - a ** 3.0)	! for 'myu_factor = 1.0e-04'

	  !a = a + tau * 3.9e+00 * a ** 2.75		! for 'myu_factor = 1.0e-03' --- LS-mode, T = 
	  !a = a + tau * 3.9e+00 * a ** 3.0			! for 'myu_factor = 1.0e-03' --- LS-mode, T = 0.19
	  !a = a + tau * 1.29e+00 * a ** 3.0			! for 'myu_factor = 1.0e-03' --- LS-mode, T = 0.19
	  !a = a + tau * 4.1e+00 * a ** 4.0			! for 'myu_factor = 1.0e-03' --- LS-mode, T = 0.12

	  !a = a + tau * 4.50e+00 * a ** 1.0			! for 'myu_factor = 1.0e-03', alpha = 1 --- HS-mode (ECMI-2016, conference presentation)
	  a = a + tau * 1.315e+00 * a ** 3.0			! for 'myu_factor = 1.0e-03', alpha = 1 --- LS-mode (ECMI-2016, conference presentation)
	  !a = a + tau * 3.00e+00 * a ** 2.0			! for 'myu_factor = 1.0e-03', alpha = 1 --- S-mode (ECMI-2016, conference presentation)
	  !a = a + tau * 1.90e+00 * a ** 2.0			! for 'myu_factor = 1.0e-03', alpha = 1 --- S-mode (ECMI-2016, conference presentation)
	  !a = a + tau * 4.0e+00 * a ** 2.0			! for 'myu_factor = 1.0e-03', alpha = 1 --- S-mode (ECMI-2016, conference presentation)
	  !a = a + tau * 1.0e+00 * (a - a ** 3.0)	! for 'myu_factor = 1.0e-04', alpha = 0 --- Temperature waves (ECMI-2016, conference presentation)
	  !a = a + tau * 4.0e+00 * (a - a ** 3.0)	! for 'myu_factor = 1.0e-04', alpha = 0 --- Temperature waves (ECMI-2016, conference presentation)
   end if

   call SwapXY()

   Ay = 0.0
   bya = 0.0

   ! Computing in "y"
   call DoInY(tau / 2.0, res)

   if (res .eq. 0) then
      pause
      stop
   end if
    
   call SwapYX()

   Ax = 0.0
   bxa = 0.0

   ! Computing in "x"
   call DoInX(tau / 2.0, res)

   if (res .eq. 0) then
      pause
      stop
   end if

   Ar = 0.0
   bra = 0.0

   ! Computing in "r"
   call DoInR(tau / 2.0, res)

!DEC$ IF DEFINED (ERROR_ENERGY)
   ! Third experiment for NMPDE (2010-2011)
   !analytic_solution = theT(1.0, 2.0, 0.0, -pi / 3.0, pi, currTime * tau)

   ! First experiment for SIMULTECH 2011
   !analytic_solution = theT(1.0, 1.0, 0.0, 0.0, pi, currTime * tau)

   ! Nonlinear problem (for SIMULTECH 2011, the third experiment; for Int. J for Numerical Methods in Engineering, the fourth experiment; NOVA Publishers, the fourth experiment; WCE2013, first experiment)
   analytic_solution = theT(9.0, 5.0 * alpha, 3.0, -2.5, 50.0, currTime * tau)
      
   ! New nonlinear problem (for Int. J for Numerical Methods in Engineering and for Springer)
   !analytic_solution = theT(7.0, 1.0, 0.5, 0.0 * pi / 4.0, 15.0, 0.7, 0.0 * pi / 4.0, 4.0, 10.0, 100.0, currTime * tau)
!DEC$ ENDIF

   ! Computing the total mass and energy
!DEC$ IF DEFINED (MPI)   
   if (processID .eq. 0) then
!DEC$ ENDIF

!DEC$ IF DEFINED (ERROR_ENERGY)
   a = a - analytic_solution
!DEC$ ENDIF
   call GetMass(mass_dummy, energy(currTime))
!DEC$ IF DEFINED (ERROR_ENERGY)
   a = a + analytic_solution
!DEC$ ENDIF
   call GetMass(mass(currTime), energy_aux)
!DEC$ IF DEFINED (ERROR_ENERGY)
   energy(currTime) = energy(currTime) / energy_aux	! Actually, this yields (numerics - analytics) / numerics, which is not exactly the relative error, as well. Indeed, it should be (numerics - analytics) / analytics instead.
!DEC$ ENDIF

   if (mod(currTime, outputStep) .eq. 0) then
      write (*, '(F15.10, A, D20.10, D20.10)') currTime * tau, ': ', mass(currTime), energy(currTime)
	  write (*, '(D20.10, A, D20.10, A, D20.10)') cflR, ' ', cflX, ' ', cflY
   end if

   if (mod(currTime, outputToDiskStep) .eq. 0) then
      cx = cx + 1
      call SaveToDisk(cx)
   end if

!DEC$ IF DEFINED (MPI)   
   end if
!DEC$ ENDIF
   
   !write (*, *) 'Press to continue ... ->'
   !pause
   !stop
end do		! End of the cycle 'currTime'

!DEC$ IF DEFINED (MPI)   
if (processID .eq. 0) then
!DEC$ ENDIF

! Saving the output
! Mass
open (1, file = 'mass.tx2')
!write (1, '(D20.10)') mass
write (1, '(E20.10E3)') mass
close (1)

! Energy
open (1, file = 'energy.tx2')
!write (1, '(D20.10)') energy
write (1, '(E20.10E3)') energy
close (1)

pause

!DEC$ IF DEFINED (MPI)   
end if
!DEC$ ENDIF

deallocate(mass, energy)

deallocate(workingArrayR, workingArrayX, workingArrayY) 
deallocate(Ay, bya, ipivy)
deallocate(Ax, bxa, ipivx)
deallocate(Ar, bra, ipivr)

!deallocate(indexSuperDiagY4, indexSuperDiagY3, indexSuperDiagY2, indexSuperDiagY, indexMainDiagY, indexSubDiagY, indexSubDiagY2, indexSubDiagY3, indexSubDiagY4)
!deallocate(indexSuperDiagX4, indexSuperDiagX3, indexSuperDiagX2, indexSuperDiagX, indexMainDiagX, indexSubDiagX, indexSubDiagX2, indexSubDiagX3, indexSubDiagX4)
deallocate(indexSuperDiagY, indexMainDiagY, indexSubDiagY)
deallocate(indexSuperDiagX, indexMainDiagX, indexSubDiagX)
deallocate(indexSuperDiagR, indexMainDiagR, indexSubDiagR)

!deallocate(a, rhs, analytic_solution)
deallocate(a)
!!!deallocate(rhs, analytic_solution)

deallocate(myu)

!DEC$ IF DEFINED (MYU_FACTOR_IS_ARRAY)
deallocate(myu_factor)
!DEC$ ENDIF

!!!deallocate(cosphi, sinphi, tanphi)

deallocate(mycosTabY, mysinTabY, mytanTabY)
deallocate(mycosTabX, mysinTabX, mytanTabX)
deallocate(my_r2)

!!!deallocate(xs, ys)

!DEC$ IF DEFINED (MPI)
call mpi_finalize(ierr)
!DEC$ ENDIF

end program Diffusion3DSphere