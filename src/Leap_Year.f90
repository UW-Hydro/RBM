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
      LOGICAL FUNCTION LEAP_YEAR (YEAR)
!
      IMPLICIT NONE
!
      INTEGER :: YEAR

      LEAP_YEAR = .FALSE.
      IF (MOD(YEAR,4) .EQ. 0)   LEAP_YEAR = .TRUE.
      IF (MOD(YEAR,100) .EQ. 0) LEAP_YEAR = .FALSE.
      IF (MOD(YEAR,400) .EQ. 0) LEAP_YEAR = .TRUE.
      RETURN
      END FUNCTION LEAP_YEAR
