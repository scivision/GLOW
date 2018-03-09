program testdrive

use cglow, only: nw, jmax,verbose,nbins

implicit none

integer,parameter :: idate=16355      
real,parameter :: ut=0,glat=80,glong=0,f107a=70,f107=70,f107p=70,ap=4,ef=1,ec=2000
real :: phi(nbins,3)

real, allocatable :: z(:),zeta(:,:)

integer :: argc,i
character(80) :: argv

argc = command_argument_count()

do i=1,argc
  call get_command_argument(i,argv)
  select case (argv)
    case ('-v','-d','--verbose','--debug')
      verbose=.true.
    case default
   
  end select
end do


jmax=102
allocate(z(jmax), zeta(nw,jmax))
! altitude grid
z = [80.,  81.,  82.,  83.,  84.,  85.,  86.,  87.,  88.,  89., &
     90.,  91.,  92.,  93.,  94.,  95.,  96.,  97.,  98.,  99., &
    100., 101., 102., 103., 104., 105., 106., 107., 108., 109., &    
    110.,111.5, 113.,114.5, 116., 118., 120., 122., 124., 126., &
    128., 130., 132., 134., 137., 140., 144., 148., 153., 158., &
    164., 170., 176., 183., 190., 197., 205., 213., 221., 229., &
    237., 245., 254., 263., 272., 281., 290., 300., 310., 320., &
    330., 340., 350., 360., 370., 380., 390., 400., 410., 420., &
    430., 440., 450., 460., 470., 480., 490., 500., 510., 520., &
    530., 540., 550., 560., 570., 580., 590., 600., 610., 620., &
    630., 640. ]

call egrid(phi(:,1),phi(:,2)) ! energy bins
call maxt (ef,ec,phi(:,1),phi(:,2),nbins,0,0,0,phi(:,3)) ! precipitation


call glowbasic(idate,ut,glat,glong,f107a,f107,f107p,ap, z,jmax,nw,nbins,phi,verbose,&
              zeta)


end program
