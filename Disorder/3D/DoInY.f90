subroutine DoInY(tau, res)

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

real, pointer :: mycosTabY(:), mysinTabY(:), mytanTabY(:)

real, pointer :: myu(:, :, :)

real R_Earth

real, pointer :: Ay(:, :), bya(:), ipivy(:)

!integer, pointer :: indexMainDiagY(:), indexSuperDiagY(:), indexSuperDiagY2(:), indexSuperDiagY3(:), indexSuperDiagY4(:), indexSubDiagY(:), indexSubDiagY2(:), indexSubDiagY3(:), indexSubDiagY4(:)
integer, pointer :: indexMainDiagY(:), indexSuperDiagY(:), indexSubDiagY(:)
real, pointer :: workingArrayY(:)

real cflR, cflX, cflY

!DEC$ IF DEFINED (MPI)
integer ierr, processID, NumberOfProcesses
!DEC$ ENDIF

!common /arhs/ a, rhs
common /a/ a
common /trigsY/ mycosTabY, mysinTabY, mytanTabY
common /steps/ dr, dx, dy
common /sizes/ lr, lx, ly
common /myu/ myu
common /Radius/ R_Earth
common /matvecY/ Ay, bya, ipivy
!common /diagsY/ indexMainDiagY, indexSuperDiagY, indexSuperDiagY2, indexSuperDiagY3, indexSuperDiagY4, indexSubDiagY, indexSubDiagY2, indexSubDiagY3, indexSubDiagY4, workingArrayY
common /diagsY/ indexMainDiagY, indexSuperDiagY, indexSubDiagY, workingArrayY
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
integer iFrom, iTo

!DEC$ IF DEFINED (PGI)
integer info
external dgetrf, dgetrs
!DEC$ ENDIF

!DEC$ IF DEFINED (MPI)
integer j
!DEC$ ENDIF

! Body
iFrom = 1
iTo = lx

!DEC$ IF DEFINED (MPI)
!!! IMPORTANT NOTE: THE NUMBER OF PROCESSES MUST BE SUCH THAT iTo - iFrom + 1 = const FOR ALL PROCESSES !!!
!!! OTHERWISE THE SUBROUTINE 'mpi_gather()' WILL PROVIDE INCORRECT RESULTS
iFrom = 1 + lx * 1.0 / NumberOfProcesses * processID
iTo = 1 + lx * 1.0 / NumberOfProcesses * (processID + 1) - 1
!write (*, *) 'DoInY: process #', processID, ': from ', iFrom, ' to ', iTo
!write (*, *) ' '
!DEC$ ENDIF

allocate(a_next(iTo - iFrom + 1, ly, lr))

do ij_r = 1, lr
   do i = iFrom, iTo
      Ay(:, :) = 0.0

      Ay(1, 1) = 1.0
      Ay(1, ly - 1) = -1.0
    
      Ay(ly, 2) = 1.0
      Ay(ly, ly) = -1.0
    
      ! Computation.
      workingArrayY = reshape(Ay, (/product(shape(Ay))/))
      !workingArrayY(indexSuperDiagY4) = 0.0
      !workingArrayY(indexSuperDiagY3) = 0.0
      !workingArrayY(indexSuperDiagY2) = 0.0
      !workingArrayY(indexSubDiagY2) = 0.0
      !workingArrayY(indexSubDiagY3) = 0.0
      !workingArrayY(indexSubDiagY4) = 0.0

!DEC$ IF DEFINED (APPROXIMATION_TYPE_1)
	  ! This code implements the wide stencil, without fictitious nodes for the diffusion coefficient
	     !workingArrayY(indexSuperDiagY2) = -myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
	     !workingArrayY(indexMainDiagY) = 1.0 / tau + (myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) + myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2)) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
	     !workingArrayY(indexSubDiagY2) = -myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)

	  ! This code implements the standard stencil, with fictitious nodes for the diffusion coefficient
	     workingArrayY(indexSuperDiagY) = -(myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) + myu(i, 2 : ly - 1, ij_r) * mycosTabY(2 : ly - 1)) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
	     workingArrayY(indexMainDiagY) = 1.0 / tau + (myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) + 2.0 * myu(i, 2 : ly - 1, ij_r) * mycosTabY(2 : ly - 1) + myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2)) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
	     workingArrayY(indexSubDiagY) = -(myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2) + myu(i, 2 : ly - 1, ij_r) * mycosTabY(2 : ly - 1)) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
!DEC$ ENDIF
!DEC$ IF DEFINED (APPROXIMATION_TYPE_2)
	     workingArrayY(indexMainDiagY) = 1.0 / tau
!DEC$ ENDIF
         Ay = reshape(workingArrayY, (/ly, ly/))
    
      if (isDiagPredom(Ay) .eq. 0) then
         write (*, *) 'The matrix Ay is not diagonally predominant.'
         !res = 0
         !return
      end if
    
	  bya(1) = -a(i, 1, ij_r) + a(i, ly - 1, ij_r)
      bya(ly) = -a(i, 2, ij_r) + a(i, ly, ij_r)
    
!DEC$ IF DEFINED (APPROXIMATION_TYPE_1)
	  ! This code implements the wide stencil, without fictitious nodes for the diffusion coefficient
	     !bya(2 : ly - 1) = a(i, 2 : ly - 1, ij_r) * (1.0 / tau - (myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) + myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2)) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)) + &
		 !										a(i, [4 : ly, 3], ij_r) * myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) + &
		 !										a(i, [ly - 2, 1 : ly - 3], ij_r) * myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2) / (8.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)

	  ! This code implements the standard stencil, with fictitious nodes for the diffusion coefficient
	     bya(2 : ly - 1) = a(i, 2 : ly - 1, ij_r) * (1.0 / tau - (myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) + 2.0 * myu(i, 2 : ly - 1, ij_r) * mycosTabY(2 : ly - 1) + myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2)) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)) + &
												a(i, 3 : ly, ij_r) * (myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) + myu(i, 2 : ly - 1, ij_r) * mycosTabY(2 : ly - 1)) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) + &
												a(i, 1 : ly - 2, ij_r) * (myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2) + myu(i, 2 : ly - 1, ij_r) * mycosTabY(2 : ly - 1)) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
!DEC$ ENDIF
!DEC$ IF DEFINED (APPROXIMATION_TYPE_2)
	  ! This code implements the wide stencil, without fictitious nodes for the diffusion coefficient
	     !bya(2 : ly - 1) = a(i, 2 : ly - 1, ij_r) * (1.0 / tau - (myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) + myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2)) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)) + &
		 !										a(i, [4 : ly, 3], ij_r) * myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) + &
		 !										a(i, [ly - 2, 1 : ly - 3], ij_r) * myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2) / (4.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)

	  ! This code implements the standard stencil, with fictitious nodes for the diffusion coefficient
	     bya(2 : ly - 1) = a(i, 2 : ly - 1, ij_r) * (1.0 / tau - (myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) + 2.0 * myu(i, 2 : ly - 1, ij_r) * mycosTabY(2 : ly - 1) + myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2)) / (2.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)) + &
												a(i, 3 : ly, ij_r) * (myu(i, 3 : ly, ij_r) * mycosTabY(3 : ly) + myu(i, 2 : ly - 1, ij_r) * mycosTabY(2 : ly - 1)) / (2.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0) + &
												a(i, 1 : ly - 2, ij_r) * (myu(i, 1 : ly - 2, ij_r) * mycosTabY(1 : ly - 2) + myu(i, 2 : ly - 1, ij_r) * mycosTabY(2 : ly - 1)) / (2.0 * mycosTabY(2 : ly - 1) * (R_Earth * dy) ** 2.0)
!DEC$ ENDIF

      cflY = max(maxval(abs(Ay)) * tau, cflY)

!DEC$ IF .NOT. DEFINED (PGI)
      a_next(i - iFrom + 1, :, ij_r) = Ay .ix. bya
!DEC$ ELSE
      ipivy(:) = 0
      call dgetrf(ly, ly, Ay, ly, ipivy, info)

      if (info .eq. 0) then
         call dgetrs('N', ly, 1, Ay, ly, ipivy, bya, ly, info)

         a_next(i - iFrom + 1, :, ij_r) = bya
      else
         write (*, *) 'The matrix Ay is singular.'
      end if
!DEC$ ENDIF

      ! Update
      a(i, :, ij_r) = a_next(i - iFrom + 1, :, ij_r)
   end do
end do

!write (*, *) 'before: ', processID, a(1, 99, 9), a(16, 99, 9), a(31, 99, 9), a(46, 99, 9)
!write (*, *) ' '

!DEC$ IF DEFINED (MPI)
do ij_r = 1, lr
   do j = 1, ly
      call mpi_gather(a_next(1, j, ij_r), iTo - iFrom + 1, MPI_REAL8, a(1, j, ij_r), iTo - iFrom + 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
   end do
end do

call mpi_bcast(a, product(shape(a)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
!DEC$ ENDIF

!write (*, *) 'after: ', processID, a(1, 99, 9), a(16, 99, 9), a(31, 99, 9), a(46, 99, 9)
!write (*, *) ' '

deallocate(a_next)

res = 1

end subroutine DoInY