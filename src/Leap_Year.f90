!*************************************************************************************
!
!  LEAP_YEAR
!
!  Author:       Dr. David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         November 1, 2002
!  
!
!  Input:
!     year  -  Gregorian year (integer)
!  Output:
!     Function return value = .true. if year is a leap year, and .false. otherwise.
!*************************************************************************************
!
      LOGICAL FUNCTION Leap_Year(nyear)
!
      IMPLICIT NONE
!
      INTEGER :: nyear
!
      leap_year = .FALSE.
      IF (MOD(nyear,4) .EQ. 0)   leap_year = .TRUE.
      IF (MOD(nyear,100) .EQ. 0) leap_year = .FALSE.
      IF (MOD(nyear,400) .EQ. 0) leap_year = .TRUE.
      RETURN
      END FUNCTION Leap_Year
