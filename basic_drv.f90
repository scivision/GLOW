program testdrive

implicit none

integer,parameter :: idate=16355        
real,parameter :: ut=0,glat=80,glong=0,f107a=70,f107=70,f107p=70,ap=4,ef=1,ec=2000

call glowbasic(idate,ut,glat,glong,f107a,f107,f107p,ap,ef,ec)


end program
