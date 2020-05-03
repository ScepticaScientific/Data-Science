subroutine DoInX(tau, res)

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

real, pointer :: mycosTabX(:), mysinTabX(:), mytanTabX(:)

real, pointer :: myu(:, :, :)

real R_Earth

real, pointer :: Ax(:, :), bxa(:), ipivx(:)

!integer, pointer :: indexMainDiagX(:), indexSuperDiagX(:), indexSuperDiagX2(:), indexSuperDiagX3(:), indexSuperDiagX4(:), indexSubDiagX(:), indexSubDiagX2(:), indexSubDiagX3(:), indexSubDiagX4(:)
integer, pointer :: indexMainDiagX(:), indexSuperDiagX(:), indexSubDiagX(:)
real, pointer :: workingArrayX(:)

real cflR, cflX, cflY

!DEC$ IF DEFINED (MPI)
integer ierr, processID, NumberOfProcesses
!DEC$ ENDIF

!common /arhs/ a, rhs
common /a/ a
common /trigsX/ mycosTabX, mysinTabX, mytanTabX
common /steps/ dr, dx, dy
common /sizes/ lr, lx, ly
common /myu/ myu
common /Radius/ R_Earth
common /matvecX/ Ax, bxa, ipivx
!common /diagsX/ indexMainDiagX, indexSuperDiagX, indexSuperDiagX2, indexSuperDiagX3, indexSuperDiagX4, indexSubDiagX, indexSubDiagX2, indexSubDiagX3, indexSubDiagX4, workingArrayX
common /diagsX/ indexMainDiagX, indexSuperDiagX, indexSubDiagX, workingArrayX
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
integer j
integer jFrom, jTo

!DEC$ IF DEFINED (PGI)
integer info
external dgetrf, dgetrs
!DEC$ ENDIF

! Body
jFrom = 1
jTo = ly

!DEC$ IF DEFINED (MPI)
!!! IMPORTANT NOTE: THE NUMBER OF PROCESSES MUST BE SUCH THAT jTo - jFrom + 1 = const FOR ALL PROCESSES !!!
!!! OTHERWISE THE SUBROUTINE 'mpi_gather()' WILL PROVIDE INCORRECT RESULTS
jFrom = 1 + ly * 1.0 / NumberOfProcesses * processID
jTo = 1 + ly * 1.0 / NumberOfProcesses * (processID + 1) - 1
!write (*, *) 'DoInX: process #', processID, ': from ', jFrom, ' to ', jTo
!write (*, *) ' '
!DEC$ ENDIF

allocate(a_next(lx, jTo - jFrom + 1, lr))

do ij_r = 1, lr
   do j = jFrom, jTo
      Ax(:, :) = 0.0

      Ax(1, 1) = 1.0
      Ax(1, lx - 1) = -1.0
    
      Ax(lx, 2) = 1.0
      Ax(lx, lx) = -1.0
    
      ! Computation.
      workingArrayX = reshape(Ax, (/product(shape(Ax))/))
      !workingArrayX(indexSuperDiagX4) = 0.0
      !workingArrayX(indexSuperDiagX3) = 0.0
      !workingArrayX(indexSuperDiagX2) = 0.0
      !workingArrayX(indexSubDiagX2) = 0.0
      !workingArrayX(indexSubDiagX3) = 0.0
      !workingArrayX(indexSubDiagX4) = 0.0

!DEC$ IF DEFINED (APPROXIMATION_TYPE_1)
	  ! This code implements the wide stencil, without fictitious nodes for the diffusion coefficient
         !workingArrayX(indexSuperDiagX2) = -myu(3 : lx, j, ij_r) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
	     !workingArrayX(indexMainDiagX) = 1.0 / tau + (myu(3 : lx, j, ij_r) + myu(1 : lx - 2, j, ij_r)) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
         !workingArrayX(indexSubDiagX2) = -myu(1 : lx - 2, j, ij_r) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)

	  ! This code implements the standard stencil, with fictitious nodes for the diffusion coefficient
         workingArrayX(indexSuperDiagX) = -(myu(3 : lx, j, ij_r) + myu(2 : lx - 1, j, ij_r)) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
	     workingArrayX(indexMainDiagX) = 1.0 / tau + (myu(3 : lx, j, ij_r) + 2.0 * myu(2 : lx - 1, j, ij_r) + myu(1 : lx - 2, j, ij_r)) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
         workingArrayX(indexSubDiagX) = -(myu(1 : lx - 2, j, ij_r) + myu(2 : lx - 1, j, ij_r)) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
!DEC$ ENDIF
!DEC$ IF DEFINED (APPROXIMATION_TYPE_2)
	     workingArrayX(indexMainDiagX) = 1.0 / tau
!DEC$ ENDIF
         Ax = reshape(workingArrayX, (/lx, lx/))
    
      if (isDiagPredom(Ax) .eq. 0) then
         write (*, *) 'The matrix Ax is not diagonally predominant.'
         !res = 0
         !return
      end if
    
      bxa(1) = -a(1, j, ij_r) + a(lx - 1, j, ij_r)		! This corresponds to the Crank-Nicolson periodic boundary conditions a^n_1 + a^{n+1}_1 = a^n_{lx-1} + a^{n+1}_{lx-1} (see the matrix elements Ax(1, 1) and Ax(1, lx - 1) above)
      bxa(lx) = -a(2, j, ij_r) + a(lx, j, ij_r)			! ...						...						...				  a^n_2 + a^{n+1}_2 = a^n_{lx} + a^{n+1}_{lx} (see the matrix elements Ax(lx, 2) and Ax(lx, lx) above)
    
!DEC$ IF DEFINED (APPROXIMATION_TYPE_1)
	  ! This code implements the wide stencil, without fictitious nodes for the diffusion coefficient
	     !bxa(2 : lx - 1) = a(2 : lx - 1, j, ij_r) * (1.0 / tau - (myu(3 : lx, j, ij_r) + myu(1 : lx - 2, j, ij_r)) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)) + &
		 ! 										a([4 : lx, 3], j, ij_r) * myu(3 : lx, j, ij_r) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) + &
		 !										a([lx - 2, 1 : lx - 3], j, ij_r) * myu(1 : lx - 2, j, ij_r) / (8.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)

      ! This code implements the standard stencil, with fictitious nodes for the diffusion coefficient
	     bxa(2 : lx - 1) = a(2 : lx - 1, j, ij_r) * (1.0 / tau - (myu(3 : lx, j, ij_r) + 2.0 * myu(2 : lx - 1, j, ij_r) + myu(1 : lx - 2, j, ij_r)) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)) + &
												a(3 : lx, j, ij_r) * (myu(3 : lx, j, ij_r) + myu(2 : lx - 1, j, ij_r)) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) + &
												a(1 : lx - 2, j, ij_r) * (myu(1 : lx - 2, j, ij_r) + myu(2 : lx - 1, j, ij_r)) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
!DEC$ ENDIF
!DEC$ IF DEFINED (APPROXIMATION_TYPE_2)
      ! This code implements the wide stencil, without fictitious nodes for the diffusion coefficient
	     !bxa(2 : lx - 1) = a(2 : lx - 1, j, ij_r) * (1.0 / tau - (myu(3 : lx, j, ij_r) + myu(1 : lx - 2, j, ij_r)) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)) + &
		 !										a([4 : lx, 3], j, ij_r) * myu(3 : lx, j, ij_r) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) + &
		 !										a([lx - 2, 1 : lx - 3], j, ij_r) * myu(1 : lx - 2, j, ij_r) / (4.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)

	  ! This code implements the standard stencil, with fictitious nodes for the diffusion coefficient
	     bxa(2 : lx - 1) = a(2 : lx - 1, j, ij_r) * (1.0 / tau - (myu(3 : lx, j, ij_r) + 2.0 * myu(2 : lx - 1, j, ij_r) + myu(1 : lx - 2, j, ij_r)) / (2.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)) + &
												a(3 : lx, j, ij_r) * (myu(3 : lx, j, ij_r) + myu(2 : lx - 1, j, ij_r)) / (2.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0) + &
												a(1 : lx - 2, j, ij_r) * (myu(1 : lx - 2, j, ij_r) + myu(2 : lx - 1, j, ij_r)) / (2.0 * (R_Earth * mycosTabX(j) * dx) ** 2.0)
!DEC$ ENDIF

      cflX = max(maxval(abs(Ax)) * tau, cflX)

!DEC$ IF .NOT. DEFINED (PGI)
      a_next(:, j - jFrom + 1, ij_r) = Ax .ix. bxa
!DEC$ ELSE
      ipivx(:) = 0
      call dgetrf(lx, lx, Ax, lx, ipivx, info)

      if (info .eq. 0) then
         call dgetrs('N', lx, 1, Ax, lx, ipivx, bxa, lx, info)

         a_next(:, j - jFrom + 1, ij_r) = bxa
      else
         write (*, *) 'The matrix Ax is singular.'
      end if
!DEC$ ENDIF
       
      ! Update
      a(:, j, ij_r) = a_next(:, j - jFrom + 1, ij_r)
   end do
end do

!write (*, *) 'before: ', processID, a(1, 1, 5), a(1, 16, 5), a(1, 31, 5), a(1, 46, 5)
!write (*, *) ' '

!DEC$ IF DEFINED (MPI)
do ij_r = 1, lr
   call mpi_gather(a_next(1, 1, ij_r), lx * (jTo - jFrom + 1), MPI_REAL8, a(1, 1, ij_r), lx * (jTo - jFrom + 1), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
end do

call mpi_bcast(a, product(shape(a)), MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
!DEC$ ENDIF

!write (*, *) 'after: ', processID, a(1, 1, 5), a(1, 16, 5), a(1, 31, 5), a(1, 46, 5)
!write (*, *) ' '

deallocate(a_next)

res = 1

end subroutine DoInX