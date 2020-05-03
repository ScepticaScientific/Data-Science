subroutine DoInR(tau, res)

!DEC$ DEFINE APPROXIMATION_TYPE_1
! APPROXIMATION_TYPE_1 - Crank-Nicolson
! APPROXIMATION_TYPE_2 - Euler (purely explicit)
! APPROXIMATION_TYPE_3 - Purely implicit

use linear_operators

use isDiagPredomFunction

implicit none

!DEC$ IF DEFINED (MPI)
include 'mpif.h'
!DEC$ ENDIF

! Variables (general)
real dr, dx, dy
integer lr, lx, ly

!real, pointer :: a(:, :, :), rhs(:, :, :)
real, pointer :: a(:, :, :)

real, pointer :: my_r2(:)

real, pointer :: myu(:, :, :)

real R_Earth

real, pointer :: Ar(:, :), bra(:), ipivr(:)

integer, pointer :: indexMainDiagR(:), indexSuperDiagR(:), indexSubDiagR(:)
real, pointer :: workingArrayR(:)

real cflR, cflX, cflY

!DEC$ IF DEFINED (MPI)
integer ierr, processID, NumberOfProcesses
!DEC$ ENDIF

!common /arhs/ a, rhs
common /a/ a
common /nontrigR/ my_r2
common /steps/ dr, dx, dy
common /sizes/ lr, lx, ly
common /myu/ myu
common /Radius/ R_Earth
common /matvecR/ Ar, bra, ipivr
common /diagsR/ indexMainDiagR, indexSuperDiagR, indexSubDiagR, workingArrayR
common /CFL/ cflR, cflX, cflY

!DEC$ IF DEFINED (MPI)
common /mpiData/ ierr, processID, NumberOfProcesses
!DEC$ ENDIF

! Variables (partial)
real tau
integer res

real, pointer :: a_next(:, :, :)
! NB: The parameter 'a_next' is actually needed only when the MPI technology is used. If so, it is needed in the subroutine 'mpi_gather()'
! as the sender buffer. If we used direcly 'a' as the sender buffer (alongside with using it as the receiver, of course), then we would obtain
! improper gathering of solutions by the process No. 0 from the others.

integer ij_r
integer i
integer j
integer jFrom, jTo

!return

! Body
jFrom = 1
jTo = ly

!DEC$ IF DEFINED (MPI)
!!! IMPORTANT NOTE: THE NUMBER OF PROCESSES MUST BE SUCH THAT jTo - jFrom + 1 = const FOR ALL PROCESSES !!!
!!! OTHERWISE THE SUBROUTINE 'mpi_gather()' WILL PROVIDE INCORRECT RESULTS
jFrom = 1 + ly * 1.0 / NumberOfProcesses * processID
jTo = 1 + ly * 1.0 / NumberOfProcesses * (processID + 1) - 1
!write (*, *) 'DoInR: process #', processID, ': from ', jFrom, ' to ', jTo
!write (*, *) ' '
!DEC$ ENDIF

allocate(a_next(lx, jTo - jFrom + 1, lr))

do i = 1, lx
   do j = jFrom, jTo
      Ar(:, :) = 0.0

      !Ar(1, 1) = 1.0
      Ar(1, 1) = -1.0
      Ar(1, 2) = 1.0
    
	  !Ar(lr, lr) = 1.0
      Ar(lr, lr - 1) = -1.0
      Ar(lr, lr) = 1.0
    
      ! Computation.
      workingArrayR = reshape(Ar, (/product(shape(Ar))/))

      !if (ApprOrder .eq. 2) then
!DEC$ IF DEFINED (APPROXIMATION_TYPE_1)
	  ! This code implements the standard stencil, with fictitious nodes for the diffusion coefficient
	     !workingArrayR(indexSuperDiagR) = -(myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 3 : lr) * my_r2(3 : lr)) / (4.0 * (R_Earth * dr) ** 2.0)
	     !workingArrayR(indexMainDiagR) = 1.0 / tau + (myu(i, j, 3 : lr) * my_r2(3 : lr) + 2.0 * myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (4.0 * (R_Earth * dr) ** 2.0)
	     !workingArrayR(indexSubDiagR) = -(myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (4.0 * (R_Earth * dr) ** 2.0)
	     workingArrayR(indexSuperDiagR) = -(myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 3 : lr) * my_r2(3 : lr)) / (4.0 * my_r2(2 : lr - 1) * dr ** 2.0)
	     workingArrayR(indexMainDiagR) = 1.0 / tau + (myu(i, j, 3 : lr) * my_r2(3 : lr) + 2.0 * myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (4.0 * my_r2(2 : lr - 1) * dr ** 2.0)
	     workingArrayR(indexSubDiagR) = -(myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (4.0 * my_r2(2 : lr - 1) * dr ** 2.0)
!DEC$ ENDIF
!DEC$ IF DEFINED (APPROXIMATION_TYPE_3)
      ! This code implements the standard stencil, with fictitious nodes for the diffusion coefficient
	     !workingArrayR(indexSuperDiagR) = -(myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 3 : lr) * my_r2(3 : lr)) / (2.0 * (R_Earth * dr) ** 2.0)
	     !workingArrayR(indexMainDiagR) = 1.0 / tau + (myu(i, j, 3 : lr) * my_r2(3 : lr) + 2.0 * myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (2.0 * (R_Earth * dr) ** 2.0)
	     !workingArrayR(indexSubDiagR) = -(myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (2.0 * (R_Earth * dr) ** 2.0)
	     workingArrayR(indexSuperDiagR) = -(myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 3 : lr) * my_r2(3 : lr)) / (2.0 * my_r2(2 : lr - 1) * dr ** 2.0)
	     workingArrayR(indexMainDiagR) = 1.0 / tau + (myu(i, j, 3 : lr) * my_r2(3 : lr) + 2.0 * myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (2.0 * my_r2(2 : lr - 1) * dr ** 2.0)
	     workingArrayR(indexSubDiagR) = -(myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (2.0 * my_r2(2 : lr - 1) * dr ** 2.0)
!DEC$ ENDIF
         Ar = reshape(workingArrayR, (/lr, lr/))
      !end if
    
      if (isDiagPredom(Ar) .eq. 0) then
         write (*, *) 'The matrix Ar is not diagonally predominant.'
         !res = 0
         !return
      end if
    
	  bra(1) = 0.0		! Zero Neumann boundary condition at the internal boundary
      bra(lr) = 0.0		! Zero Neumann boundary condition at the external boundary
    
      !if (ApprOrder .eq. 2) then
!DEC$ IF DEFINED (APPROXIMATION_TYPE_1)
      ! This code implements the standard stencil, with fictitious nodes for the diffusion coefficient
	     !bra(2 : lr - 1) = a(i, j, 2 : lr - 1) * (1.0 / tau - (myu(i, j, 3 : lr) * my_r2(3 : lr) + 2.0 * myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (4.0 * (R_Earth * dr) ** 2.0)) + &
		 !									a(i, j, 3 : lr) * (myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 3 : lr) * my_r2(3 : lr)) / (4.0 * (R_Earth * dr) ** 2.0) + &
		 !									a(i, j, 1 : lr - 2) * (myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (4.0 * (R_Earth * dr) ** 2.0)
	     bra(2 : lr - 1) = a(i, j, 2 : lr - 1) * (1.0 / tau - (myu(i, j, 3 : lr) * my_r2(3 : lr) + 2.0 * myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (4.0 * my_r2(2 : lr - 1) * dr ** 2.0)) + &
	 	  									a(i, j, 3 : lr) * (myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 3 : lr) * my_r2(3 : lr)) / (4.0 * my_r2(2 : lr - 1) * dr ** 2.0) + &
	 	 									a(i, j, 1 : lr - 2) * (myu(i, j, 2 : lr - 1) * my_r2(2 : lr - 1) + myu(i, j, 1 : lr - 2) * my_r2(1 : lr - 2)) / (4.0 * my_r2(2 : lr - 1) * dr ** 2.0)
!DEC$ ENDIF
!DEC$ IF DEFINED (APPROXIMATION_TYPE_3)
      ! This code implements the standard stencil, with fictitious nodes for the diffusion coefficient
	     bra(2 : lr - 1) = a(i, j, 2 : lr - 1) / tau
!DEC$ ENDIF
      !end if

      cflR = max(maxval(abs(Ar)) * tau, cflR)

      a_next(i, j - jFrom + 1, :)  = Ar .ix. bra

      ! Update
      a(i, j, :) = a_next(i, j - jFrom + 1, :)
   end do
end do

!write (*, *) 'before: ', processID, a(1, 1, 1), a(1, 31, 1)
!write (*, *) ' '

!DEC$ IF DEFINED (MPI)
do ij_r = 1, lr
   call mpi_gather(a_next(1, 1, ij_r), lx * (jTo - jFrom + 1), MPI_REAL8, a(1, 1, ij_r), lx * (jTo - jFrom + 1), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
end do

call mpi_bcast(a, product(shape(a)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
!DEC$ ENDIF

!write (*, *) 'after: ', processID, a(1, 1, 1), a(1, 31, 1)
!write (*, *) ' '

deallocate(a_next)

res = 1

end subroutine DoInR