! Subroutine EGRID sets up electron energy grid

! This software is part of the GLOW model.  Use is governed by the open source
! academic research license agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon, 1/1992
! Refactored to f90, SCS, 6/2016

!   nbins   number of bins in the electron energy grid
! Outputs:
!   ener    energy at center of each bin, eV
!   del     width of each bin, eV

 
subroutine egrid (ener, del)
  use cglow,only: nbins

  implicit none

  real,intent(out) :: ener(nbins), del(nbins)

  integer :: n

  do n=1,nbins
    if (n <= 21) then
      ener(n) = 0.5 * float(n)
    else
      ener(n) = exp (0.05 * float(n+26))
    endif
  enddo

  del(1) = 0.5
  del(2:) = ener(2:)-ener(:nbins-1)

  ener = ener - del / 2.0

end subroutine egrid
