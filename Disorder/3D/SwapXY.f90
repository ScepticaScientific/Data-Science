subroutine SwapXY()

implicit none

! Variables (general)
!real, pointer :: a(:, :, :), rhs(:, :, :)
real, pointer :: a(:, :, :)

!!!real, pointer :: xy_r(:), x(:), y(:)

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
integer i

! Body
llxNew = pi
llyNew = 2.0 * pi

!allocate(xNew(int((llxNew - dx) / dx + 1)));					lxNew = size(xNew, 1)
lxNew = int((llxNew - dx) / dx + 1)
!!!xNew = [0 : lxNew - 1] * dx + dx / 2.0

!allocate(yNew(int((llyNew + dy) / dy + 1)));					lyNew = size(yNew, 1)
lyNew = int((llyNew + dy) / dy + 1)
!!!yNew = -lly + dy / 2.0 + [0 : lyNew - 1] * dy

allocate(myuNew(lxNew, lyNew, lr))
allocate(aNew(lxNew, lyNew, lr))

do ij_r = 1, lr
   do i = 1, lxNew
      myuNew(i, :, ij_r) = [myu(i, :, ij_r), myu(size(myu, 1) - 1 - i, size(myu, 2) : 1 : -1, ij_r), myu(i, 1 : 2, ij_r)]
    
      aNew(i, :, ij_r) = [a(i, :, ij_r), a(size(a, 1) - 1 - i, size(a, 2) : 1 : -1, ij_r), a(i, 1 : 2, ij_r)]
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
!do i = 1, lxNew
!   un(i, :) = unNew(i, :)
!   vn(i, :) = vnNew(i, :)
!   zn(i, :) = znNew(i, :)
!end do
! The assignement "un = unNew" (and the similars) does not work when "lxNew" and "lyNew" are large (say, for lxNew = 360 and lyNew = 180)
! So, we used a cycle

deallocate(a)
allocate(a(lxNew, lyNew, lr))
a = aNew
!do i = 1, lxNew
!   U(i, :) = UNew(i, :)
!   V(i, :) = VNew(i, :)
!   a(i, :) = aNew(i, :)
!   zb(i, :) = zbNew(i, :)
!end do

!!!deallocate(xNew)
!!!deallocate(yNew)
deallocate(myuNew)
deallocate(aNew)

end subroutine SwapXY
