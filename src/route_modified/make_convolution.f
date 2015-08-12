
      SUBROUTINE MAKE_CONVOLUTION
     & (Nheat,NCOL, NROW, NOB, PMAX, DAYS, CATCHIJ 
     & ,BASE, RUNO, FLOW, Qin,KE, UH_DAY, UH_S, FRACTION
     & ,XC, YC, SIZE, DPREC, FLOWPATH,GRID_CELL,clen,ICOL
     & ,NDAY,IDAY,IMONTH,IYEAR,MO, YR, NYR,Ntest)

      IMPLICIT NONE

      INTEGER     N, I, J, DAYS, NDAY, II, JJ, Nheat  
      INTEGER     NCOL,NROW,ICOL,NOB,PMAX,KE,UH_DAY
      INTEGER     CATCHIJ(PMAX,2)
      INTEGER     NYR,iwrite,Ntest
      REAL        UH_S(PMAX,KE+UH_DAY-1)
c
c     Modified by JRY 9/30/2009 to include Qin from VIC_Cell
c
      REAL        BASE(DAYS), RUNO(DAYS), FLOW(DAYS),Qin(days) 
      REAL        FRACTION(NCOL,NROW)

      REAL        PI, RERD, FACTOR, FACTOR_SUM

      LOGICAL TORF

      PARAMETER   (RERD  = 6371229.0)    !radius of earth in meters
c 
      CHARACTER*80 GRID_CELL,FLOWPATH

      INTEGER DPREC, CLEN

      REAL        JLOC, ILOC
      REAL        XC, YC, SIZE
      REAL        AREA, AREA_SUM

      REAL        STORAGE, K_CONST
      REAL        DUM1(9),DUM2,DUM3,DUM4,DUM5
      REAL        DUM6,DUM7,DUM8,DUM9,DUM10
  
      INTEGER     ii1,IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
      INTEGER MO(12*NYR),YR(12*NYR)
      CHARACTER*20 LOC   ! Yixin hacked

C     *** 0 <= K_CONST = 1.0 
C *** K_CONST smaller 1.0 makes it a simple linear storage

      K_CONST = 1.0

      PI = ATAN(1.0) * 4.0

      AREA_SUM   = 0.0
      FACTOR_SUM = 0.0
      DO I = 1,NDAY
         FLOW(I) = 0.0
      END DO

      DO N = 1,NOB       !is this the gridcell loop?
         STORAGE = 0.0
         DO I = 1,NDAY
            RUNO(I) = 0.0
            BASE(I) = 0.0
         END DO
         II = CATCHIJ(N,1)
         JJ = CATCHIJ(N,2)
c
c     the grid has been flipped left to right
c     find the revised cooordinates

         ILOC=XC + (ICOL-II)*SIZE + SIZE/2.0
         JLOC=YC + JJ*SIZE - SIZE/2.0

C        CONVERSIONFACTOR for mm/day to ft**3/sec         <--?????

         AREA =  RERD**2*ABS(SIZE)*PI/180*             !give area of box in 
     &        ABS(SIN((JLOC-SIZE/2.0)*PI/180)-         !square meters
     $        SIN((JLOC+SIZE/2.0)*PI/180))
         
         AREA_SUM = AREA_SUM + AREA ! AREA_SUM: not used - Yixin
         iwrite=icol-ii+1
c
         FACTOR = FRACTION(II,JJ)*35.315*AREA/(86400.0*1000.0)  !convert to sq.mi. - NOT CORRECT! - Yixin
                                                                !&mult. by cell fract
                                                                !Convert flow
                                                                !from [mm/day] to
                                                                ! [cfs] - Yixin
c 
         FACTOR_SUM = FACTOR_SUM + FACTOR  ! FACTOR_SUM: not used - Yixin
c 
! ============= Yixin hacked ===================
         call create_vic_names(jloc,iloc,loc,clen,dprec)


   !      INQUIRE(FILE=TRIM(FLOWPATH)//TRIM(GRID_CELL),
   !  $                     EXIST=TORF)
         INQUIRE(FILE=TRIM(FLOWPATH)//LOC(1:CLEN),
     $                     EXIST=TORF)
   !    if (.not.TORF) write(91,*) 'No file '
   !  &  ,TRIM(FLOWPATH)//TRIM(GRID_CELL)
       if (.not.TORF) write(91,*) 'No file '
     &  ,TRIM(FLOWPATH)//LOC(1:CLEN)
c 
         if(torf)then
c
   !        OPEN(20,FILE=TRIM(FLOWPATH)//TRIM(GRID_CELL)
   !  $         ,STATUS='OLD',ERR=9001)
           OPEN(20,FILE=TRIM(FLOWPATH)//LOC(1:CLEN)
     $         ,STATUS='OLD',ERR=9001)
! ============= Yixin hacked (end) ===================
c
c       Read the header on the VIC fluxes files
c 
           do i=1,6
              read(20,*)
           end do 
c
c     read VIC model output: <year> <month> <day> <p> <et> <runoff> <baseflow>
c
           DO I = 1,NDAY
             READ(20,*,END=9001,ERR=9001) IYEAR(I),IMONTH(I),IDAY(I)
     &                                   ,DUM2, DUM3, RUNO(I), BASE(I)
       if (runo(i).lt.0.0001.and.base(i).lt.0.0001) write(94,*)
c
c     check to be sure dates in VIC file start at same time specified
c     in input file
             if(I.eq.1) then
               if(IYEAR(I).ne.YR(1) . or. IMONTH(I).ne.MO(1)) then
                  print*, 'VIC output file does not match specified '
                  print*, 'period in input file.'
                  stop
               endif
             endif
           END DO
	 else
C
	   do i=1,nday
             IYEAR(I)=9999
             IMONTH(I)=99
             IDAY(I)=99
	     runo(i)=0
	     base(i)=0
	   end do
	 endif

         DO I = 1,NDAY
            RUNO(I) = RUNO(I) * FACTOR
            BASE(I) = BASE(I) * FACTOR
         END DO
         DO I = 1,NDAY
c
c  Add Qin  JRY 9/30/2009 
c    
         Qin(i) = 0.0
            DO J = 1,KE+UH_DAY-1
               IF ((I-J+1) .GE. 1) THEN
c
c  Added by JRY 9/30/2009
c
c                 Qin(i) = Qin(i)+UH_S(N,J)*(BASE(I-J+1)+RUNO(I-J+1))
                 FLOW(I) = FLOW(I)+UH_S(N,J)*(BASE(I-J+1)+RUNO(I-J+1))
                END IF
            END DO
c            
c  Added by JRY 9/30/2009             
c            
c            FLOW(I) = FLOW(I)+Qin(i)
         END DO
         CLOSE(20)
      END DO
      RETURN
 9001 WRITE(95,*) 'Error reading time-series data, ',
     $     'insufficient data or missing input file',
     $     TRIM(FLOWPATH)//TRIM(GRID_CELL)
      END
