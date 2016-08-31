
      SUBROUTINE READ_DIREC(DIREC,NCOL,NROW,H,XC,YC,SIZE
     $     ,FILENAME,IROW,ICOL)

c  reads the flow direction file.

      IMPLICIT NONE

      INTEGER NCOL,NROW,I,J,IROW,ICOL,IMISS
      INTEGER DIREC(NCOL,NROW,2) 
      INTEGER H(NCOL,NROW)
      REAL XC, YC, SIZE
      CHARACTER*72 FILENAME
      CHARACTER*14 CDUM 

      OPEN(10, FILE = FILENAME, STATUS='OLD',ERR=9001)

      READ(10,*) CDUM, ICOL    !note: ARC/INFO style header
      READ(10,*) CDUM, IROW
      READ(10,*) CDUM, XC
      READ(10,*) CDUM, YC
      READ(10,*) CDUM, SIZE
      READ(10,*) CDUM, IMISS
      print*, 'Flow direction file dimensions:'
      print*, 'cols:', icol,' rows:',irow
      print*, 'x-corner:', xc,' ycorner:',yc
      print*, 'size:', size, ' void:', imiss
      IF(IROW.GT.NROW .OR. ICOL.GT.NCOL)THEN
         WRITE(*,*) 'Incorrect dimensions:'
         WRITE(*,*) 'Reset nrow and ncol in main to;',
     $        irow, icol
         STOP
      ENDIF
      
      DO J = IROW,1,-1
         READ(10,*) (H(I,J), I=ICOL,1,-1) 
      END DO      
      CLOSE(10)

      DO I = 1, ICOL
         DO J = 1,IROW
            IF (H(I,J) .EQ. 9) THEN     !note outlet is 9, not 0
               DIREC(I,J,1) = 0
               DIREC(I,J,2) = 0
            ELSE IF (H(I,J) .EQ. 1) THEN
               DIREC(I,J,1) = I
               DIREC(I,J,2) = J+1
            ELSE IF (H(I,J) .EQ. 2) THEN
               DIREC(I,J,1) = I-1
               DIREC(I,J,2) = J+1
            ELSE IF (H(I,J) .EQ. 3) THEN
               DIREC(I,J,1) = I-1
               DIREC(I,J,2) = J
            ELSE IF (H(I,J) .EQ. 4) THEN 
               DIREC(I,J,1) = I-1
               DIREC(I,J,2) = J-1
            ELSE IF (H(I,J) .EQ. 5) THEN
               DIREC(I,J,1) = I
               DIREC(I,J,2) = J-1
            ELSE IF (H(I,J) .EQ. 6) THEN
               DIREC(I,J,1) = I+1
               DIREC(I,J,2) = J-1
            ELSE IF (H(I,J) .EQ. 7) THEN
               DIREC(I,J,1) = I+1
               DIREC(I,J,2) = J
            ELSE IF (H(I,J) .EQ. 8) THEN
               DIREC(I,J,1) = I+1
               DIREC(I,J,2) = J+1
            END IF
         END DO
      END DO
      RETURN
 9001 WRITE(*,*) 'CANNOT OPEN INPUT FILE IN READ_DIREC',
     $  FILENAME
      STOP
      END

