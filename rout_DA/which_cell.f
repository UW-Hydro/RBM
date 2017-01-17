	SUBROUTINE WHICH_CELL(XC,YC,SIZE,IROW,ICOL,LAT,LON,PI,PJ)
c  search the grid to find the cell corresponding to the input latitude
c  and longitude

c     Corrected type of XRC to real and corrected the position of XRC 
c     and YLC so that lat lon increment for I and J is correct.
c     A. Hamlet 4/27/01
  
      real LAT,LON,SIZE,XC,YC,XRC
      integer IROW,ICOL,PI,PJ

c     calculate the cell indexes (PI,PJ) for the latitude and longitude
      
      XRC = XC + SIZE*ICOL + SIZE/2
      PI = -1
      PJ = -1
      DO I = 1, ICOL  
        DO J = 1, IROW
          IF(I .eq. ICOL .and. J .eq. IROW) THEN
            IF(LON .le. XRC-(I-1)*SIZE .and. LON .ge. XRC-I*SIZE .and.
     &       LAT .ge. (J-1)*SIZE+(YC-SIZE/2) .and. LAT .le. J*SIZE+
     &  (YC-SIZE/2)) THEN	  
              PI = I
              PJ = J
            END IF
          ELSE IF(I .lt. ICOL .and. J .eq. IROW) THEN
            IF(LON .le. XRC-(I-1)*SIZE .and. LON .gt. XRC-I*SIZE .and.
     &       LAT .ge. (J-1)*SIZE+(YC-SIZE/2) .and. LAT .le. J*SIZE+
     &  (YC-SIZE/2)) THEN
              PI = I
              PJ = J
            END IF
          ELSE
            IF(LON .le. XRC-(I-1)*SIZE .and. LON .gt. XRC-I*SIZE .and.
     &       LAT .ge. (J-1)*SIZE+(YC-SIZE/2) .and. LAT .lt. J*SIZE+
     & (YC-SIZE/2)) THEN

c	       print*, XRC,I,J
c	       print*, XRC-(I-1)*SIZE, XRC-I*SIZE,(J-1)*SIZE+YC


              PI = I
              PJ = J
            END IF
          END IF
        END DO
      END DO
      IF (PI .eq. -1) THEN
        PRINT*, 'NO SUBBASIN OUTLET CELL FOUND IN GRID'
        PRINT*, 'CHECK YOUR LAT-LON VALUES IN THE INPUT FILE'
        STOP
      END IF

      RETURN
      END
