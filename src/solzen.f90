module sun_angles
  use cglow, only: wp,pi
  implicit none

contains
! Subroutine SOLZEN

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file glowlicense.txt.
! For more information see the file glow.txt.

! Stan Solomon, 1988.
! Temporary Y2K fix-up, SCS, 2005.
! Refactored to f90, SCC, 2016.

! Returns Solar Zenith Angle SZA in degrees for specified date in form yyddd or yyyyddd,
! universal time in seconds, geographic latitude and longitude in degrees.

subroutine solzen (idate, ut, glat, glong, sza)

  integer,intent(in) :: idate
  real(wp),intent(in) :: ut, glat, glong
  real(wp),intent(out) :: sza

  real(wp) :: rlat, rlong, sdec, srasn, gst, rh, cossza

  rlat = glat * pi/180._wp
  rlong = glong * pi/180._wp
  call suncor (idate, ut, sdec, srasn, gst)
  rh = srasn - (gst+rlong)
  cossza = sin(sdec)*sin(rlat) + cos(sdec)*cos(rlat)*cos(rh)
  sza = acos(cossza) * 180._wp/pi

end subroutine solzen


! Subroutine SUNCOR returns the declination SDEC and right ascension
! SRASN of the sun in GEI coordinates, radians, for a given date IDATE
! in yyddd or yyyyddd format, universal time UT in seconds, and Greenwich sidereal
! time GST in radians.  Reference:  C.T. Russell, Geophysical Coordinate Transforms.

subroutine suncor (idate, ut, sdec, srasn, gst)

  integer,intent(in) :: idate
  real(wp),intent(in) :: ut
  real(wp),intent(out) :: sdec, srasn, gst

  real(wp) :: fday, dj, t, vl, g, slong, obliq, slp, sind, cosd
  integer :: iyr, iday


  fday=ut/86400._wp
  iyr=idate/1000
  iday=idate-iyr*1000

! Temporary Y2K fix-up:
! Should work with either yyddd or yyyyddd format from 1950 to 2050.
! Note deteriorating accuracy after ~2050 anyway.
! Won't work after 2100 due to lack of a leap year.

  if (iyr >= 1900) iyr=iyr-1900
  if (iyr < 50) iyr=iyr+100

  dj=365*iyr+(iyr-1)/4+iday+fday-0.5_wp
  t=dj/36525._wp
  vl = mod(279.696678_wp+.9856473354_wp*dj,360._wp)
  gst = mod(279.696678_wp+.9856473354_wp*dj+360._wp*fday+180._wp,360._wp) * pi/180._wp
  g = mod(358.475845_wp+.985600267_wp*dj,360._wp) * pi/180._wp
  slong=vl+(1.91946_wp-.004789_wp*t)*sin(g)+.020094_wp*sin(2._wp*g)
  obliq=(23.45229_wp-0.0130125_wp*t) *pi/180._wp
  slp=(slong-.005686_wp) * pi/180._wp
  sind=sin(obliq)*sin(slp)
  cosd=sqrt(1._wp-sind**2._wp)
  sdec=atan(sind/cosd)
  srasn=pi-atan2(1._wp/tan(obliq)*sind/cosd, -cos(slp)/cosd)

end subroutine suncor

end module sun_angles
