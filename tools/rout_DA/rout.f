
      PROGRAM rout
c
c     Routing algorithm developed by D. Lohmann.
c
c     Modified to allow more flexible array dimensions and
c     the removal of harcoded file names.
c
c     Code maintained by G. O'Donnell (tempgd@hydro.washington.edu)
c     See WA Hydrology Homepage for operational details.

c     Modified 5/99 to read in the uh_s array if it has already
c     been generated in a previous run. 

c     Modified 2/2001 by edm to include month and year in output
c     and also check dates in VIC output files and calculate NDAYS
c
      IMPLICIT NONE
      integer IARGC

      integer isaleap,julian
      external isaleap,julian
 
c     change dimensions here
c     nrow and ncol should be larger than the grid
c     nyr should equal run length yrs+1
      INTEGER NROW, NCOL, DAYS, NYR, Nflow, Nheat
      PARAMETER (NROW = 195, NCOL = 245)
      PARAMETER (NYR = 100)

c     no changes after here
      REAL    DT
      INTEGER KE, LE, TMAX, UH_DAY, PMAX
      PARAMETER (DAYS=NYR*366)
      PARAMETER (KE   = 12    )
      PARAMETER (LE   = 48    )
      PARAMETER (DT   = 3600.0)
      PARAMETER (UH_DAY = 96  )
      PARAMETER (TMAX = UH_DAY*24)
      PARAMETER (PMAX = 5000   )
      
      INTEGER DIREC(NCOL,NROW,2)
      REAL    VELO(NCOL,NROW), DIFF(NCOL,NROW) 
      REAL    XMASK(NCOL,NROW), FRACTION(NCOL,NROW)
      REAL    UH_BOX(PMAX,KE), UHM(NCOL,NROW,LE)
      REAL    UH_S(PMAX,KE+UH_DAY-1)
      real    a_d(ncol,nrow),b_d(ncol,nrow)
     $       ,a_w(ncol,nrow),b_w(ncol,nrow)
     $       ,aa_d,bb_d,aa_w,bb_w
     $       ,a_ddum,b_ddum,a_wdum,b_wdum
c
c
c  Add Qin - JRY- 9/30/2009
c
      REAL    BASE(DAYS), RUNO(DAYS), FLOW(DAYS), Qin(days)

      INTEGER IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
      INTEGER MO(12*NYR),YR(12*NYR)
      INTEGER NO_OF_BOX
      INTEGER CATCHIJ(PMAX,2)
      INTEGER H(NCOL,NROW)
      
      INTEGER PI, PJ
      REAL    UH_DAILY(PMAX,UH_DAY)
      REAL    FR(TMAX,2)
c
c     ndd,Nout,nrec,nseg added by JRY - 10/30/2009
c      
      INTEGER Flow_Cells,Force_Cells,NR,ndd,Nout,nrec,nseg
      INTEGER IROW, ICOL
      INTEGER LP,M,Y
      INTEGER J

      CHARACTER*80 UH_STRING,UH_DRCTRY         ! new, AW
      CHARACTER*25 loc,NAME5,NAME         !was 21
c      CHARACTER*24  NAME5       !was  5 6/25/2009
      CHARACTER*72 FILE_INPUT, FILENAME
      character*10 lake10
c
c     Deleted by JRY -11/4/2009
c
c      CHARACTER*72 INPATH, OUTPATH
      CHARACTER*72 FLOWPATH, HEATPATH
      CHARACTER*72 FLOWOUT, HEATOUT

      INTEGER DPREC,clen,fdleni,fdleno,hdleni,hdleno
     &       ,Nflow_pd, Nheat_pd,Ntest
      REAL    AREA, FACTOR_SUM

      REAL XC, YC, SIZE

      REAL FDUM
      LOGICAL TORF,rbm10,write_flow
      INTEGER NDAY
      INTEGER NMONTHS

C**** variables for monthly means *****************************

      INTEGER DAYS_IN_MONTH(12)
      DATA DAYS_IN_MONTH /31,28,31,30,31,30,31,31,30,31,30,31/

      INTEGER START_YEAR, STOP_YEAR, FIRST_YEAR, LAST_YEAR
      INTEGER START_MO, STOP_MO, FIRST_MO, LAST_MO
      integer start_day,stop_day
      integer strtmet_yr,strtmet_mo,strtmet_dy
     &       ,stopmet_yr,stopmet_mo,stopmet_dy
     &       ,skip_MET
      integer MET_strt,MET_stop,FLO_strt,FLO_stop
      REAL MONTHLY(12*NYR)
      REAL MONTHLY_mm(12*NYR)
      REAL YEARLY(12)
      REAL YEARLY_mm(12)

C***********************************************************
C     OPEN NECESSARY FILES
C***********************************************************
      fraction(ncol,nrow)=9999.
      fraction(ncol-1,nrow-1)=8888.
c      write(*,*) 'Ntest'
c      read(*,*) Ntest
c     INPUT FILE READ
c     process commandline args
      IF(IARGC().NE.1)THEN
         PRINT*, 'USAGE:  rout <infile>'
         STOP
      ENDIF
      CALL GETARG(1,FILE_INPUT)
      OPEN(1,FILE=FILE_INPUT,STATUS='OLD',ERR=9001)
      READ(1,'(//A)') FILENAME
      CALL READ_DIREC(DIREC,NCOL,NROW,H,XC,YC,SIZE,
     $     FILENAME,IROW,ICOL)
c     velo file
      READ(1,*)
      READ(1,*) TORF
      IF(TORF)THEN
         READ(1,'(A)') FILENAME
         CALL READ_VELO(VELO,NCOL,NROW,FILENAME,
     $        IROW,ICOL)
      ELSE
         READ(1,*) FDUM
         CALL INIT_ARRAY(VELO,NCOL,NROW,FDUM)
      ENDIF
c     diff file
      READ(1,*)
      READ(1,*)TORF
      IF(TORF)THEN
         READ(1,'(A)') FILENAME
         CALL READ_DIFF(DIFF,NCOL,NROW,FILENAME,
     $        IROW,ICOL)
      ELSE
         READ(1,*) FDUM
         CALL INIT_ARRAY(DIFF,NCOL,NROW,FDUM)
      ENDIF
c     xmask file
      READ(1,*)
      READ(1,*)TORF
      IF(TORF)THEN
         READ(1,'(A)') FILENAME
         CALL READ_XMASK(XMASK,NCOL,NROW,FILENAME,
     $        IROW,ICOL)
      ELSE
         READ(1,*) FDUM
         CALL INIT_ARRAY(XMASK,NCOL,NROW,FDUM)
      ENDIF
c     read fraction file
      write(*,*) 'Reading fraction file'
      READ(1,*)
      READ(1,*)TORF
      IF(TORF)THEN
         READ(1,'(A)') FILENAME
         CALL READ_FRACTION(FRACTION,NCOL,NROW,FILENAME,
     $        IROW,ICOL)
      ELSE
         READ(1,*) FDUM
         CALL INIT_ARRAY(FRACTION,NCOL,NROW,FDUM)
      ENDIF
c
c     read file with Leopold relationships
c     (UW_JRY_2010/12/14)
c
      READ(1,*) lake10
      write(*,*) 'FDUM = ',lake10
      READ(1,*)TORF
      IF(TORF)THEN
c
c   Reading the Leopold coefficients from a file (e.g. Impound.Out
c   Note to Michelle/Wietse: This can be modified without too much
c   effort to read each of the four coefficients from a separate
c   file. If you do, you must also modify the subroutine 
c   READ_LEOPOLD in the read_routines.f file.   (JRY 3/15/2011)
c
         READ(1,'(A)') FILENAME
         CALL READ_LEOPOLD(a_d,b_d,a_w,b_w,NCOL,NROW,FILENAME,
     $        IROW,ICOL)
      ELSE
         READ(1,*) a_ddum,b_ddum,a_wdum,b_wdum
         CALL INIT_ARRAY(a_d,NCOL,NROW,a_ddum)
         CALL INIT_ARRAY(b_d,NCOL,NROW,b_ddum)
         CALL INIT_ARRAY(a_w,NCOL,NROW,a_wdum)
         CALL INIT_ARRAY(b_w,NCOL,NROW,b_wdum)
      ENDIF
c
c     station file
c
      READ(1,'(/A)')FILENAME
      write(*,*) 'Station file',FILENAME
      OPEN(10,FILE=FILENAME)
c     read input path and precision of VIC filenames
c      READ(1,'(/A)')INPATH
c
c     If RBM10 = .TRUE. open forcing directory and set up
c     direct access files
c
      read(1,*) rbm10
      write(*,*) 'RBM10?',rbm10
c     Input files for flow
c
      READ(1,'(/A)') FLOWPATH
      write(*,*) 'Flowpath',FLOWPATH
      fdleni=index(FLOWPATH,' ')-1
      write(*,*) 'fdleni',fdleni
c
c     Input files for heat flux
c
      if (rbm10) READ(1,'(A)')HEATPATH
      hdleni=index(HEATPATH,' ')-1
      write(*,*) 'hdleni ',hdleni
c
c     Read precision of file extension
c
      READ(1,*)DPREC,Nflow_pd,Nheat_pd
c
c     output pathname for flows
c
c      READ(1,'(/A)')OUTPATH
      READ(1,'(/A)') FLOWOUT
      fdleno=index(FLOWOUT,' ')-1
      if (rbm10) then
      write(*,*) 'DA files '
     & ,flowout(1:fdleno)
c
c     Direct access file opened for writing results used by RBM10 - JRY 10/30/2009
c

        open (15,FILE=FLOWOUT(1:fdleno),FORM='FORMATTED'
     &          ,ACCESS='DIRECT',RECL=60)
c
c    Output pathname for heat flux
c
        READ(1,'(A)') HEATOUT
        hdleno=index(HEATOUT,' ')-1
        write(*,*) 'DA files ' 
     &             ,heatout(1:hdleno)                  
c        
      open(16,FILE=HEATOUT(1:hdleno),FORM='FORMATTED'
     &          ,ACCESS='DIRECT',RECL=50)
      end if      
        
c
c     number of days to process
      READ(1,*)
c     start and end year/month from VIC simulation
      READ(1,*) START_YEAR, START_MO, start_day
     &         ,STOP_YEAR, STOP_MO, stop_day
c
c     Calculate Julian date of start and stop for VIC flow
c

      FLO_strt=julian(start_year,start_mo,start_day)
      FLO_stop=julian(stop_year,stop_mo,stop_day)
      nday=FLO_stop-FLO_strt+1
      write(*,*) 'VIC flow Julian,#days -',FLO_strt,nday
  
c     calculate number of days & months in simulation
      M=START_MO
      Y=START_YEAR
      NMONTHS = 0
      NDAY=0
      DO J=START_MO,12*(STOP_YEAR-START_YEAR)+STOP_MO
        IF(M.EQ.2) THEN
           LP=isaleap(Y)
        ELSE 
           LP=0
        ENDIF
        NDAY = NDAY+DAYS_IN_MONTH(M)+LP
        NMONTHS = NMONTHS + 1
        MO(NMONTHS) = M
        YR(NMONTHS) = Y
        M = M + 1
        IF (M .GT. 12) THEN
            M = 1
            Y  = Y + 1
        ENDIF
      END DO
      IF(NDAY.GT.DAYS) THEN
         PRINT*, 'IN ROUT.F RESET DAYS TO ', NDAY
         STOP
      ENDIF
      PRINT*,'Simulation Days NDAY = ',NDAY, ' NMONTHS = ',NMONTHS
c
c     start and end day for full_data record
c
      read(1,*) strtmet_yr,strtmet_mo,strtmet_dy
     &         ,stopmet_yr,stopmet_mo,stopmet_dy
c
c     Julian date of beginning and end of meteorology record (full_data)
c
      MET_strt=julian(strtmet_yr,strtmet_mo,strtmet_dy)
      MET_stop=julian(stopmet_yr,stopmet_mo,stopmet_dy)
      write(*,*) 'MET_strt',MET_strt,'MET_stop',MET_stop
c
c     Numberof Met records to skip
c 
      skip_MET=Nheat_pd*(FLO_strt-MET_strt)
      write(*,*) 'Number of MET records to skip',skip_MET
c
c     start and end year/month for writing output
c 
      READ(1,*) FIRST_YEAR, FIRST_MO, LAST_YEAR, LAST_MO
c
c     Read Unit Hydrograph
c 
      READ(1,'(/A)')FILENAME
c
c     Read directory to store final *uh_s files
c 
      READ(1,'(/A)') UH_DRCTRY

C***********************************************************

      CALL MAKE_UHM(UHM,VELO,DIFF,XMASK,NCOL,NROW,LE,DT,
     $        IROW,ICOL)

C     Loop over required stations
      Nflow=0
      Nheat=0
      read(10,*) Flow_Cells,Force_Cells
 100  CONTINUE
      READ(10,*,END=110) 
     &     NR, NSEG, NAME, PI, PJ, AREA
      clen=INDEX(NAME,' ')-1
      Nheat=Nheat+1
      write_flow=.false.
      IF (NR .EQ. 1) THEN
         WRITE(*,'(2I2,2X,A,I4,I4,G12.6)') 
     &        NR, NSEG, NAME, PI, PJ
      READ(10,'(A80)',END=110) UH_STRING   !new, AW:  uh_string
c
c     note, the arrays are flipped left to right
c
c
c     Set Leopold Coefficients (UW_JRY_2011/02/03)
c 
        aa_d=a_d(PI,PJ)
        bb_d=b_d(PI,PJ)
        aa_w=a_w(PI,PJ)
        bb_w=b_w(PI,PJ)
        if (bb_w.lt.0.01)
     &  write(70,*) ' Rout ',PJ,PI,aa_d,bb_d,aa_w,bb_w
c 
         Nflow=Nflow+1
         write_flow=.true.
         PI=ICOL+1-PI
         NAME5 = NAME
         loc=name
c 
         CALL SEARCH_CATCHMENT
     &        (PI,PJ,DIREC,NCOL,NROW,
     &        NO_OF_BOX,CATCHIJ,PMAX,IROW,ICOL)
c 
         CALL READ_GRID_UH
     &        (UH_BOX, KE, PMAX, NO_OF_BOX, CATCHIJ,FILENAME)
c 
         CALL MAKE_GRID_UH
     &        (DIREC, NO_OF_BOX, UH_DAY, TMAX, PI, PJ, LE, UH_DAILY, KE,
     &        CATCHIJ,UHM, FR, PMAX, NCOL, NROW, UH_BOX, UH_S,
     &        UH_STRING,UH_DRCTRY,NAME5,clen)        !new, AW:  added uh_string
c 
         CALL MAKE_CONVOLUTION
     & (NCOL, NROW, NO_OF_BOX, PMAX, DAYS, CATCHIJ, 
     &  BASE, RUNO, FLOW, KE, UH_DAY, UH_S, FRACTION, FACTOR_SUM,
     &  XC, YC, SIZE, DPREC, FLOWPATH,ICOL,NDAY,IDAY,IMONTH,IYEAR,
     &  MO, YR, NYR)
c 
         end if
c 
         print*, 'writing data...',nheat
         CALL WRITE_DATA
c
c  Add Qin - JRY- 9/30/2009
c 
     &     (write_flow,rbm10,FLOW, Qin, DAYS,NAME5
     &     ,flowout,fdleno,heatpath,hdleni,LOC,clen
     &     ,Flow_Cells,Force_Cells,Nflow,Nflow_pd,Nheat,Nheat_pd
     &     ,nday,IDAY,IMONTH,IYEAR,skip_MET
     &     ,aa_d,bb_d,aa_w,bb_w)

c
c

      GOTO 100
 110  CONTINUE
      
      STOP
 9001 WRITE(*,*) 'CANNOT OPEN: ', FILE_INPUT 
      END
c     ***********************************************
c     FUNCTION  ISALEAP

      integer function isaleap( iyr )

c     return 1 if a leap yr else 0

      if( (mod(iyr,4) .eq. 0 .and. mod(iyr,100) .ne.0)
     $                       .or. mod(iyr,400) .eq. 0) then
         isaleap = 1
      else
         isaleap = 0
      endif

      end
c
c
      INTEGER FUNCTION Julian (YEAR,MONTH,DAY)
C
C---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
C   DATE (YEAR,MONTH,DAY).
C
      INTEGER YEAR,MONTH,DAY,I,J,K
C
      I= YEAR
      J= MONTH
      K= DAY
C
      Julian=
     1   K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12)
     2  /12-3*((I+4900+(J-14)/12)/100)/4
C
      RETURN
      END

