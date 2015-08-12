      SUBROUTINE MAKE_CONVOLUTION
     & (NCOL, NROW, NOB, PMAX, CATCHIJ, 
     &  FRACTION, FACTOR_SUM,
     $     XC, YC, SIZE, DPREC, INPATH,ICOL,AREA_ARG)

      IMPLICIT NONE

      INTEGER     N, II, JJ  
      INTEGER     NCOL,NROW,ICOL,NOB,PMAX
      INTEGER     CATCHIJ(PMAX,2)
      REAL        FRACTION(NCOL,NROW)

      REAL        PI, RERD, FACTOR, FACTOR_SUM,AREA_ARG

      PARAMETER   (RERD  = 6371229.0)    !radius of earth in meters

      CHARACTER*20 LOC
      CHARACTER*72 INPATH

      INTEGER DPREC, CLEN

      REAL        JLOC, ILOC
      REAL        XC, YC, SIZE
      REAL        AREA, AREA_SUM

      PI = ATAN(1.0) * 4.0

      AREA_SUM   = 0.0
      FACTOR_SUM = 0.0

      DO N = 1,NOB       !is this the gridcell loop?

         II = CATCHIJ(N,1)
         JJ = CATCHIJ(N,2)
         
c     the grid has been flipped left to right
c     find the revised cooordinates

         ILOC=XC + (ICOL-II)*SIZE + SIZE/2.0
         JLOC=YC + JJ*SIZE - SIZE/2.0

C        CONVERSIONFACTOR for mm/day to ft**3/sec         <--?????


         AREA =  RERD**2*ABS(SIZE)*PI/180*             !give area of box in 
     &        ABS(SIN((JLOC-SIZE/2.0)*PI/180)-         !square meters
     $        SIN((JLOC+SIZE/2.0)*PI/180))

         AREA_SUM = AREA_SUM + AREA

c         WRITE(*,*) N, ILOC, JLOC

         FACTOR = FRACTION(II,JJ)*35.315*AREA/(86400.0*1000.0)  !convert to sq.mi.
                                                                !&mult. by cell fract
         FACTOR_SUM = FACTOR_SUM + FACTOR

c         print*, 'cellfract: ',fraction(ii,jj),' sq.mi.: ',factor,
c     &     ' tot: ',factor_sum

c         CALL CREATE_VIC_NAMES(jloc,iloc,loc,clen,dprec)

c        print*, INPATH(1:INDEX(INPATH,' ')-1)//LOC(1:CLEN)

      END DO
      AREA_ARG = factor_sum

      RETURN
      END

