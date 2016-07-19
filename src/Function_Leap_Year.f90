Program Test_Leap
implicit none
integer    :: year
logical    :: leap
10 continue
read(*,*) year
if (leap(year)) then
  write(*,*) year
  go to 10
end if
stop
end program Test_leap
!************************************************************************************************
!  LEAP
!
!  Author:       Dr. David G. Simpson
!                NASA Goddard Space Flight Center
!                Greenbelt, Maryland  20771
!
!  Date:         November 1, 2002
!
!  Input:
!     year  -  Gregorian year (integer)
!  Output:
!     Function return value = .true. if year is a leap year, and .false. otherwise.
!************************************************************************************************

      FUNCTION LEAP (YEAR) RESULT (LEAPFLAG)

      IMPLICIT NONE

      INTEGER :: YEAR
      LOGICAL :: LEAPFLAG

      LEAPFLAG = .FALSE.
      IF (MOD(YEAR,4) .EQ. 0)   LEAPFLAG = .TRUE.
      IF (MOD(YEAR,100) .EQ. 0) LEAPFLAG = .FALSE.
      IF (MOD(YEAR,400) .EQ. 0) LEAPFLAG = .TRUE.
      RETURN
      END FUNCTION LEAP
