! Subroutine ROUT writes model atmosphere and excitation rates to an output file
! in order to transfer them to Randy Gladstone's REDISTER radiative transfer program.

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Scott Baily and Stan Solomon, 9/1994
! Replaced 834 with LBH, SCS, 2/2003
! Reduced cascade contribution to 1356, SCS, 9/2003
! Included radiative recombination in 1356, SCS, 9/2003
! Refactored to f90, SCS, 12/2016
! Changed 1356, 1304, and LBH to use volume emission rate arrays, SCS, 12/2016


SUBROUTINE ROUT(ROFILE,LUN,EF,EZ,ITAIL,FRACO,FRACO2,FRACN2)

  use cglow,only: jmax,idate
  use cglow,only: ut,glat,glong,f107,f107p,f107a
  use cglow,only: zz,aglw,sza,dip,xuvfac,ztn,zti,zte,zo,zo2,zns,zn2,ecalc,zxden,zeta

  implicit none

  integer,intent(in) :: lun,itail
  real,intent(in) :: ef,ez,fraco,fraco2,fracn2
  character(*),intent(in) :: rofile

  real, dimension(jmax) :: z, zhe, e1356, e1304, e1027, e989, elbh
  integer :: j

    z(:)=zz(:)/1.e5
    zhe(:)=0.
    e1356(:)=zeta(13,:)
    e1304(:)=zeta(14,:)
    e1027(:)=aglw(7,1,:)
    e989(:)=aglw(8,1,:)
    elbh(:)=zeta(12,:)

  open(unit=lun,file=rofile,status='unknown')

  write(lun,"('   JMAX ','   SZA  ','   UT   ','   IDATE','   LAT  ','   LONG ','    DIP ')")
  write(lun,"(i8,f8.2,f8.1,i8,3f8.2)") jmax,sza*180./3.14159,ut,idate,glat,glong,dip
  write(lun,"('  F107  ','  F107p ','  F107a ',' XUVfac ')")
  write(lun,"(4f8.2)") f107,f107p,f107a,xuvfac
  write(lun,"(' Eflux  ',' Ezero  ',' Itail  ',' FracO  ',' FracO2  ',' FracN2  ')")
  write(lun,"(f8.2,f8.1,i8,3f8.2)") ef, ez, itail, fraco, fraco2, fracn2
  write(lun,"(' Alt    Tn    Ti    Te      O        O2       N2     &
  He       N        Ne       O+      1356     1304     1027      989&
  LBH')")

  do j=1,jmax
    write(lun,"(0p,f6.1,3f6.0,1p,12e9.2)") &
      z(j),ztn(j),zti(j),zte(j),zo(j),zo2(j),zn2(j),zhe(j), &
      zns(j),ecalc(j),zxden(3,j),e1356(j),e1304(j),e1027(j),e989(j),elbh(j)
  enddo

  close(lun)

end subroutine rout
