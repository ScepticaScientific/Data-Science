subroutine SwapYX()

implicit none

! Variables (general)
!real, pointer :: a(:, :, :), rhs(:, :, :)
real, pointer :: a(:, :, :)

real, pointer :: xy_r(:), x(:), y(:)

real dr, dx, dy
integer lr, lx, ly
real llr, llx, lly							! These are Lx and Ly, respectively, in MatLab.

real, pointer :: myu(:, :, :)

real pi

!common /arhs/ a, rhs
common /a/ a
!!!common /axes/ xy_r, x, y
common /steps/ dr, dx, dy
common /limits/ llr, llx, lly
common /sizes/ lr, lx, ly
common /myu/ myu
common /consts/ pi

! Variables (partial)
!!!real, pointer :: xNew(:)
integer lxNew
real llxNew

!!!real, pointer :: yNew(:)
integer lyNew
real llyNew							! This is LyNew in MatLab

real, pointer :: myuNew(:, :, :)
real, pointer :: aNew(:, :, :)

integer ij_r
integer j

! Body
llxNew = 2.0 * pi
llyNew = pi / 2.0
!llxNew = 1.0
!llyNew = 0.25

!allocate(xNew(int((llxNew + dx) / dx + 1)));					lxNew = size(xNew, 1)
lxNew = int((llxNew + dx) / dx + 1)
!!!xNew = [0 : lxNew - 1] * dx + dx / 2.0
	
!allocate(yNew(int((2.0 * llyNew - dy) / dy + 1)));				lyNew = size(yNew, 1)
lyNew = int((2.0 * llyNew - dy) / dy + 1)
!!!yNew = -llyNew + dy / 2.0 + [0 : lyNew - 1] * dy

allocate(myuNew(lxNew, lyNew, lr))
allocate(aNew(lxNew, lyNew, lr))

do ij_r = 1, lr
   do j = 1, (lxNew - 2) / 2
      myuNew(j, :, ij_r) = myu(j, 1 : (size(myu, 2) - 2) / 2, ij_r)
      myuNew(size(myuNew, 1) - 1 - j, :, ij_r) = myu(j, size(myu, 2) - 2 : (size(myu, 2) - 2) / 2 + 1 : -1, ij_r)
      myuNew(size(myuNew, 1) - 1 : size(myuNew, 1), :, ij_r) = myuNew(1 : 2, :, ij_r)

      aNew(j, :, ij_r) = a(j, 1 : (size(a, 2) - 2) / 2, ij_r)
      aNew(size(aNew, 1) - 1 - j, :, ij_r) = a(j, size(a, 2) - 2 : (size(a, 2) - 2) / 2 + 1 : -1, ij_r)
      aNew(size(aNew, 1) - 1 : size(aNew, 1), :, ij_r) = aNew(1 : 2, :, ij_r)
   end do
end do

llx = llxNew
lly = llyNew

!!!deallocate(x)
!!!allocate(x(lxNew))
!!!x = xNew
lx = lxNew

!!!deallocate(y)
!!!allocate(y(lyNew))
!!!y = yNew
ly = lyNew

deallocate(myu)
allocate(myu(lxNew, lyNew, lr))
myu = myuNew

deallocate(a)
allocate(a(lxNew, lyNew, lr))
a = aNew

!!!deallocate(xNew)
!!!deallocate(yNew)
deallocate(myuNew)
deallocate(aNew)

end subroutine SwapYX
