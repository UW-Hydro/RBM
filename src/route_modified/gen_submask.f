      PROGRAM gen_submask
c
c  This program takes your flow direction file, your maskfile, and an input file 
c  in which you specify a subbasin outlet latitude and longitude, and spits out 
c  a maskfile for that basin.
c
c  written by A. Wood mostly from bits of D. Lohmann and G. O'Donnell's routing code

      IMPLICIT NONE

      integer IARGC
c     change dimensions here
c     nrow and ncol should be larger than the grid
      INTEGER NROW, NCOL, CNT,FLAG
      PARAMETER (NROW = 200, NCOL = 200)
      CHARACTER CHARDUM*14
      REAL VOIDVAL,FACTOR_SUM,XC,YC,SIZE,LAT,LON

      INTEGER PMAX,I,J,L,IROW,ICOL,PI,PJ,NO_OF_BOX,DPREC

      PARAMETER (PMAX = 4000)  !Put something larger than your basin cellno here
      
      INTEGER DIREC(NCOL,NROW,2)
      REAL    XMASK(NCOL,NROW), FRACTION(NCOL,NROW)
      REAL    NEWMASK(NCOL,NROW),AREA
      INTEGER CATCHIJ(PMAX,2), H(NCOL,NROW)
      
      CHARACTER*21 NAME
      CHARACTER HEADER(5)*40
      CHARACTER*72 FILE_INPUT, FILENAME, INPATH

C***********************************************************
C     OPEN NECESSARY FILES
C***********************************************************

c     process commandline args
      IF(IARGC().NE.1)THEN
         PRINT*, 'USAGE:  gen_submask <infile>'
         PRINT*, 'infile has following format:'
	 PRINT*, '----------------------------'
         PRINT*, '# comment'
         PRINT*, '# comment'
         PRINT*, 'filename (with path) of flowdirection file'
         PRINT*, '# comment'
         PRINT*, 'filename (with path) of mask file'
         PRINT*, '# comment'
         PRINT*, 'filename (with path) of flow fraction file'
         PRINT*, '# comment'
         PRINT*, 'basinname latitude longitude (in decimal degrees)'
         PRINT*, '# comment'
         PRINT*, 
     & 'inputfile path, w/ root of numerical suffixes,ie,"fluxes_"'
         PRINT*, 'precision'
         PRINT*, '# comment'
         PRINT*, 'name of output subbasin mask file'
         STOP
      ENDIF
      CALL GETARG(1,FILE_INPUT)
      OPEN(1,FILE=FILE_INPUT,STATUS='OLD',ERR=9001)

c     flow direction file
      READ(1,'(//A)') FILENAME
      CALL READ_DIREC(DIREC,NCOL,NROW,H,XC,YC,SIZE,
     $     FILENAME,IROW,ICOL)

c     xmask file
      READ(1,'(/A)') FILENAME
      CALL READ_XMASK(XMASK,NCOL,NROW,FILENAME,
     $        IROW,ICOL,HEADER,CHARDUM,VOIDVAL)
c    passing back header information along with grid

c     read fraction file
      READ(1,'(/A)') FILENAME
      CALL READ_FRACTION(FRACTION,NCOL,NROW,FILENAME,
     $        IROW,ICOL)

c     subbasin outlet information
      READ(1,*)
      READ(1,*) NAME, LAT,LON
      PRINT*, 'submask for: ',NAME,' lat:',LAT,' lon',LON

c     read input path and precision of VIC filenames
      READ(1,'(/A)')INPATH
      READ(1,*)DPREC

c     read new xmask filename
      READ(1,'(/A)') FILENAME

C***********************************************************

c     note, the arrays are flipped left to right
c      PI=ICOL+1-PI

c     find the grid cell corresponding to the lat long
      CALL WHICH_CELL(XC,YC,SIZE,IROW,ICOL,LAT,LON,PI,PJ)

      CALL SEARCH_CATCHMENT
     &   (PI,PJ,DIREC,NCOL,NROW,
     &   NO_OF_BOX,CATCHIJ,PMAX,IROW,ICOL,NEWMASK,VOIDVAL)
         
      CALL MAKE_CONVOLUTION
     &   (NCOL, NROW, NO_OF_BOX, PMAX, 
     &   CATCHIJ, FRACTION,
     &   FACTOR_SUM,XC,YC,SIZE,DPREC,INPATH,ICOL,AREA)

c     write new subbasin mask in both mask and column format
c     mask style
      OPEN(20,file=FILENAME,status='unknown')      
      DO L = 1,5
        WRITE(20,'(A40)') HEADER(L)
      END DO
      WRITE(20,'(A14,F6.0)') CHARDUM,VOIDVAL
      DO J = IROW,1,-1
        WRITE(20,'(500F9.1)') (NEWMASK(I,J), I=ICOL,1,-1)
      END DO
      CLOSE(20)
c     column style
c      OPEN(120,file="col"//FILENAME,status='unknown')      
c      CNT = 0
c      DO J = IROW,1,-1
c        DO I = ICOL,1,-1
c          IF(XMASK(I,J) .NE. VOIDVAL) THEN
c            CNT = CNT+1
c            IF(NEWMASK(I,J) .NE. VOIDVAL) THEN
c              FLAG = 1
c              WRITE(120,'(I2,I6)') FLAG,CNT
c            ELSE
c              FLAG = 0
c              WRITE(120,'(I2,I6)') FLAG,CNT
c            END IF
c          END IF
c        END DO
c      END DO
c      CLOSE(120)   

      PRINT*, '-----------------------------------------------------'
      PRINT*, 'Subbasin outlet in row',PJ,' from bottom'
      PRINT*, 'and column',ICOL+1-PI,' from left (west)'
      PRINT*, ' ' 
      PRINT*, 'Your estimated basin area (sq. mi.) is: ', AREA

      STOP
 9001 WRITE(*,*) 'CANNOT OPEN: ', FILE_INPUT 
      END
