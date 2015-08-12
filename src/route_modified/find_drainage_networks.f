
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
      PARAMETER (NROW = 200, NCOL = 200)
      CHARACTER CHARDUM*14
      REAL VOIDVAL,FACTOR_SUM,XC,YC,SIZE,LAT,LON
      REAL CLAT,CLON
      REAL ESTAREA,CHECK
      INTEGER FPI,FPJ,FNCELLS,ST(200,200)
      INTEGER PMAX,I,J,L,IROW,ICOL,PI,PJ,NO_OF_BOX,DPREC

      PARAMETER (PMAX = 10000)  !Put something larger than your basin cellno here
      
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
      OPEN(1, FILE='input.psound',STATUS='OLD', ERR=9001)


c     get file name for output row and column list 
      READ(1,*) FILENAME
      OPEN(2,FILE=FILENAME,STATUS='UNKNOWN',ERR=9002)
c      OPEN(3,FILE='stations.cc.old',STATUS='UNKNOWN',ERR=9002)
c      OPEN(4,FILE='num_cells_upstream',STATUS='UNKNOWN', ERR=9002)


c     get data from flow direction file
      READ(1,*) FILENAME
      CALL READ_DIREC(DIREC,NCOL,NROW,H,XC,YC,SIZE,
     $     FILENAME,IROW,ICOL)


c     read fraction file
      READ(1,*) FILENAME
      CALL READ_FRACTION(FRACTION,NCOL,NROW,FILENAME,
     $        IROW,ICOL)


c     read number of locations from input file
      READ(1,*) NLOCS


c     loop over number of locations to process

      DO K=1,NLOCS

      READ(1,*) CROW, CCOL

c     search catchment expects row and column in reverse order
c     i.e. position 1,1  is IROW, ICOL
      PJ = IROW + 1 -CROW
      PI = ICOL +1 -CCOL

c     find cells routed to the cell
c     note that this version of search catchment writes a list of 
c     cells to a hard coded file

      CALL SEARCH_CATCHMENT
     &   (PI,PJ,DIREC,NCOL,NROW,
     &   NO_OF_BOX,CATCHIJ,PMAX,IROW,ICOL,NEWMASK,VOIDVAL)

      END DO


      STOP
 9001 WRITE(*,*) 'CANNOT OPEN INPUT FILE: '
 9002 WRITE(*,*) 'CANNOT OPEN OUTPUT FILE ' 
      END
