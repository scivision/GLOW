program glowbasic

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Version 0.981, 6/2017

! Adapted from glowdriver by Stan Solomon, 2/2016

! Basic single-processor driver for the GLOW model.
! Uses MSIS/IRI for input.
! Runs GLOW for designated inputs once, or multiple times.
! MPI and netCDF libraries not required.

! For definitions of use-associated variables, see subroutine GLOW and module CGLOW.

! Other definitions:
! f107p   Solar 10.7 cm flux for previous day
! ap      Ap index of geomagnetic activity
! z       altitude array, km

! Array dimensions:
! jmax    number of altitude levels
! nbins   number of energetic electron energy bins
! lmax    number of wavelength intervals for solar flux
! nmaj    number of major species
! nst     number of states produced by photoionization/dissociation
! nei     number of states produced by electron impact
! nex     number of ionized/excited species
! nw      number of airglow emission wavelengths
! nc      number of component production terms for each emission

  use, intrinsic:: iso_fortran_env, only: error_unit

  use cglow,only: jmax,nbins,lmax,nmaj,nei,nex,nw,nc,nst, idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec, &
    iscale,jlocal,kchem,xuvfac,  zz,zo,zn2,zo2,zns,znd,zno,ztn,ze,zti,zte, &
    ener,del,phitop,tir,&
    ecalc,zxden,zeta, cglow_init, data_dir

  implicit none

  character(*), parameter :: iri90_dir = 'data/iri90/'

  real,allocatable :: z(:)                    ! glow height coordinate in km (jmax)
  real,allocatable :: zun(:), zvn(:)          ! neutral wind components (not in use)
  real,allocatable :: pedcond(:), hallcond(:) ! Pederson and Hall conductivities in S/m (mho)
  real,allocatable :: outf(:,:)               ! iri output (11,jmax)
  real :: stl,fmono,emono
  integer :: j,ii,itail
  integer :: instance,iostatus
  
  data_dir    = 'data/'


! Initialize standard switches:

  iscale=1
  xuvfac=3.
  kchem=4
  jlocal=0
  itail=0
  fmono=0.
  emono=0.

! Set number of altitude levels:
  jmax = 102

! Allocate local arrays:
  allocate(z(jmax))
  allocate(zun(jmax))
  allocate(zvn(jmax))
  allocate(pedcond(jmax))
  allocate(hallcond(jmax))
  allocate(outf(11,jmax))

! Call CGLOW_INIT (module CGLOW) to set array dimensions and allocate use-associated variables:
! (This was formerly done using common blocks, including common block /cglow/.)
!
  call cglow_init
!
! Call EGRID to set up electron energy grid:
!
  call egrid (ener, del)
!
! Loop to call GLOW for designated inputs until end-of-file or any character on standard input:
!
  do instance=1,10000
!
! Get input values:
!
    write(error_unit,"('Enter date, UT, lat, lon, F107a, F107, F107p, Ap, Ef, Ec')")
    read(5,*,iostat=iostatus) idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
    if (iostatus /= 0) stop
!
! Calculate local solar time:
!
    stl = ut/3600. + glong/15.
    if (stl < 0.) stl = stl + 24.
    if (stl >= 24.) stl = stl - 24.
!
! Call MZGRID to use MSIS/NOEM/IRI inputs on default altitude grid:
!
    call mzgrid (jmax,nex,idate,ut,glat,glong,stl,f107a,f107,f107p,ap,iri90_dir, &
                 z,zo,zo2,zn2,zns,znd,zno,ztn,zun,zvn,ze,zti,zte,zxden)
!
! Call MAXT to put auroral electron flux specified by namelist input into phitop array:
!
    phitop(:) = 0.
    if (ef>.001 .and. ec>1.) call maxt (ef,ec,ener,del,nbins,itail,fmono,emono,phitop)
!
! Fill altitude array, converting to cm:
!
    zz(:) = z(:) * 1.e5     ! km to cm at all jmax levels
!
! Call GLOW to calculate ionized and excited species, airglow emission rates,
! and vertical column brightnesses:
!
    call glow
!
! Call CONDUCT to calculate Pederson and Hall conductivities:
!
    do j=1,jmax
      call conduct (glat, glong, z(j), zo(j), zo2(j), zn2(j), &
                    zxden(3,j), zxden(6,j), zxden(7,j), ztn(j), zti(j), zte(j), &
                    pedcond(j), hallcond(j))
    enddo
!
! Output section:
!
    write(error_unit,"(1x,i7,9f8.1)") idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec
    write(error_unit,"('   Z     Tn       O        N2        NO      Ne(in)      "//&
             "Ne(out)  Ionrate      O+       O2+      NO+       N(2D)    Pederson   Hall')")
    write(error_unit,"(1x,0p,f5.1,f6.0,1p,12e10.2)") (z(j),ztn(j),zo(j),zn2(j),zno(j),ze(j), &
      ecalc(j),tir(j),zxden(3,j),zxden(6,j),zxden(7,j),zxden(10,j),pedcond(j),hallcond(j),j=1,jmax)
    write(error_unit,"('   Z      3371    4278    5200    5577    6300    7320   "//&
             "10400    3644    7774    8446    3726    LBH     1356    1493    1304')")
    write(error_unit,"(1x,f5.1,15f8.2)")(z(j),(zeta(ii,j),ii=1,15),j=1,jmax)

  enddo

end program glowbasic
