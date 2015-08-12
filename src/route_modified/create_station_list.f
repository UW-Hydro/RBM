
      PROGRAM create_station_list
c
c  This program takes a flow direction file, a maskfile, and an input file 
c  containing paths and a list of station names, lat lons, and basin areas and writes
c  a station list for the routing program.  The program searches an adjustable range of cells
c  around the closest lat lon match and selects the closest basin area match.
c
c  The program will also write out a file with number of cells upstream for checking the 
c  locations. 

c  Note that the simple search only identifies an approximate location, which must be checked.
c
c  written by A. Wood mostly from bits of D. Lohmann and G. O'Donnell's routing code
c  modified by A. Hamlet 4/26/01


      IMPLICIT NONE

      
c     change dimensions here
c     nrow and ncol should be larger than the grid
      INTEGER NROW, NCOL, CNT,FLAG,NLOCS,Z,K,W,M,N
      PARAMETER (NROW = 200, NCOL = 250)
      CHARACTER CHARDUM*14
      REAL VOIDVAL,FACTOR_SUM,XC,YC,SIZE,LAT,LON
      REAL CLAT,CLON
      REAL ESTAREA,CHECK
      INTEGER FPI,FPJ,FNCELLS,ST(200,250)
      INTEGER PMAX,I,J,L,IROW,ICOL,PI,PJ,NO_OF_BOX,DPREC

      PARAMETER (PMAX = 50000)  !Put something larger than your basin cellno here
      
      INTEGER DIREC(NCOL,NROW,2)
      REAL    XMASK(NCOL,NROW), FRACTION(NCOL,NROW)
      REAL    NEWMASK(NCOL,NROW),AREA
      INTEGER CATCHIJ(PMAX,2), H(NCOL,NROW)
      
      CHARACTER*21 NAME
      CHARACTER HEADER(5)*40
      CHARACTER*72 FILE_INPUT, FILENAME, INPATH
      CHARACTER*21 NAME2,NAME3

C***********************************************************
C     OPEN NECESSARY FILES
C***********************************************************
      W=-99
      Z=1


c     open input file
      OPEN(1, FILE='input',STATUS='OLD', ERR=9001)


c     get file name for output station list and open file
      READ(1,*) FILENAME
      OPEN(2,FILE=FILENAME,STATUS='UNKNOWN',ERR=9002)
      OPEN(3,FILE='stations.cc.old',STATUS='UNKNOWN',ERR=9002)
      OPEN(4,FILE='num_cells_upstream',STATUS='UNKNOWN', ERR=9002)
      OPEN(20,FILE='Headwaters',STATUS='UNKNOWN', ERR=9002)


c     get data from flow direction file
      READ(1,*) FILENAME
      CALL READ_DIREC(DIREC,NCOL,NROW,H,XC,YC,SIZE,
     $     FILENAME,IROW,ICOL)


c     read fraction file
      READ(1,*) FILENAME
      CALL READ_FRACTION(FRACTION,NCOL,NROW,FILENAME,
     $        IROW,ICOL)


c     #############################################################
c     get the number of cells upstream of each location in the mask
c     and print out from NW to SE.  Note: change write statement below.
      
c      DO PJ=IROW,1,-1
c       DO PI=ICOL,1,-1
      
c      no_of_box=0
c     find number of cells routed to the cell
c      if(h(PI,PJ).gt.0) CALL SEARCH_CATCHMENT
c     &   (PI,PJ,DIREC,NCOL,NROW,
c     &   NO_OF_BOX,CATCHIJ,PMAX,IROW,ICOL,NEWMASK,VOIDVAL)

c      ST(PJ,PI)=NO_OF_BOX
c      if(no_of_box.eq.1) write(20,*) PJ,PI
c  100 continue
c      print*, PI,PJ,ST(PJ,PI)
      
c      END DO
c      END DO

c     print out number of cells upstream to file NW to SE


c      DO PJ=IROW,1,-1
         

c         WRITE(4,'(245I6)') (ST(PJ,PI), PI=ICOL,1,-1)

c         END DO

c     #############################################################


c     read number of locations to write to output file
      READ(1,*) NLOCS


c     loop over number of locations to process

      DO K=1,NLOCS


c     get subbasin outlet information from input file
      READ(1,*) NAME, CLAT,CLON, ESTAREA
      PRINT*, 'submask for: ',NAME,' lat:',CLAT,' lon',CLON, ESTAREA


c     note, the arrays are flipped left to right
c     i.e. column position from left=ICOL+1-PI

c      print*, XC,YC,SIZE,IROW,ICOL,CLAT,CLON 


c     search over X grid cells centered on the original lat lon values
c     choosing the one with the number of upstream cells closest to the
c     estimated basin area supplied in the input file. The search radius should probably be
c     adjusted to the resolution of the mask file.

      CHECK=1E6

      DO M=-2,2,1
         DO N= -2,2,1

            LAT=CLAT + M*SIZE
            LON=CLON + N*SIZE



c     find the grid cell corresponding to the lat long
      CALL WHICH_CELL(XC,YC,SIZE,IROW,ICOL,LAT,LON,PI,PJ)
      
      print*, 'col=',PI,'row=',PJ

      VOIDVAL = -9999.0

c     find number of cells routed to the cell
      CALL SEARCH_CATCHMENT
     &   (PI,PJ,DIREC,NCOL,NROW,
     &   NO_OF_BOX,CATCHIJ,PMAX,IROW,ICOL,NEWMASK,VOIDVAL)

c     get basin area

      CALL MAKE_CONVOLUTION
     &   (NCOL, NROW, NO_OF_BOX, PMAX, 
     &   CATCHIJ, FRACTION,
     &   FACTOR_SUM,XC,YC,SIZE,DPREC,INPATH,ICOL,AREA)

      print*, 'area = ',AREA 
           
c     store lowest difference from cell estimate in input file

      IF(ABS(ESTAREA-AREA) .lt. CHECK) THEN
         CHECK=ABS(ESTAREA-AREA)
         FPI=PI
         FPJ=PJ
         FNCELLS=NO_OF_BOX
         END IF

      
      END DO
      END DO



c     write info to station file

      NAME2 = NAME(1:5)//'.uh_s'
      NAME3 = 'NONE'

      WRITE(2,'(I1,1X,A5,1X,I3,1X,I3,1X,I4)') Z,NAME,ICOL+1-FPI,FPJ,W
      WRITE(2,'(A10)') NAME2

      WRITE(3,'(I1,1X,A5,1X,I3,1X,I3,1X,I4)') Z,NAME,ICOL+1-FPI,FPJ,W
      WRITE(3,'(A10)') NAME3
      


c     write check to standard output
      PRINT*, NAME, ICOL+1-FPI,FPJ, FNCELLS

      END DO

      STOP
 9001 WRITE(*,*) 'CANNOT OPEN INPUT FILE: '
 9002 WRITE(*,*) 'CANNOT OPEN OUTPUT FILE ' 
      END
