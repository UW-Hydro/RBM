
      SUBROUTINE SEARCH_CATCHMENT
     & (PI,PJ,DIREC,NCOL,NROW,NO_OF_BOX,CATCHIJ,PMAX,
     $  IROW,ICOL,NEWMASK,VOIDVAL)

      IMPLICIT NONE

      INTEGER PI,PJ,I,J,NCOL,NROW,PMAX,ICOL,IROW,N
      INTEGER II, JJ, III, JJJ,NO_OF_BOX
      INTEGER DIREC(NCOL,NROW,2)
      INTEGER CATCHIJ(PMAX,2)
      REAL NEWMASK(NCOL,NROW),VOIDVAL

C****** CATCHMENTS ***************************************

      NO_OF_BOX = 0

      DO I = 1, ICOL
         DO J = 1, IROW
            N=0
            NEWMASK(I,J) = VOIDVAL 
c            print*, 'i=',i,'j=',j,'no_of_box=', no_of_box
            II = I
            JJ = J


 300        CONTINUE
            N= N+1

            IF(N.GT.7000) THEN
               print*, 'sink encountered'
               print*, 'offending cell', II,JJ
               GOTO 310              
               END IF

c           print*, II,JJ, DIREC(II,JJ,1), DIREC(II,JJ,2)

            IF ((II .GT. ICOL) .OR. (II .LT.1) .OR. 
     &          (JJ .GT. IROW) .OR. (JJ .LT.1)) THEN
               GOTO 310
            END IF
            IF ((II .EQ. PI) .AND. (JJ .EQ. PJ)) THEN 
               NO_OF_BOX = NO_OF_BOX + 1
              CATCHIJ(NO_OF_BOX,1) = I
               CATCHIJ(NO_OF_BOX,2) = J
               NEWMASK(I,J) = 1.0                 !put current cell in
               GOTO 310                           !the mask
            ELSE IF ((DIREC(II,JJ,1).NE.0) .AND.    !check if the current
     &             (DIREC(II,JJ,2) .NE.0)) THEN   !ii,jj cell routes down
                     III = DIREC(II,JJ,1)         !to the subbasin outlet
                     JJJ = DIREC(II,JJ,2)         !point, following the

            IF(DIREC(III,JJJ,1).EQ.II .AND. DIREC(III,JJJ,2).EQ.JJ) THEN
               print*, 'cell ', 129-II,97-JJ,' directs to itself'
              GOTO 310
               END IF

                     II  = III                    !direction of direc(,)
                     JJ  = JJJ                    !from each cell
                     GOTO 300
             END IF                               !if you get there,
c                                                 !no_of_box increments  
 310        CONTINUE                              !and you try another   
         END DO                                   !cell.
      END DO
 
c      WRITE(*,*) 'Number of grid cells upstream of present station',
c     $     no_of_box

      RETURN
      END




