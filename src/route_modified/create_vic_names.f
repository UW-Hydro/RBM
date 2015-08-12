      SUBROUTINE CREATE_VIC_NAMES( JLOC, ILOC, EXTEN, CLEN, DPREC )


c     create string containing vic file names to be
c     appended to path given in input file

c     filenames allowed a maximum of 5 decimal places

      IMPLICIT NONE

      CHARACTER*10 JICHAR(2)
      CHARACTER*20 EXTEN
      REAL JLOC, ILOC
      INTEGER NSPACE, CLEN, CLEN_OLD, DPREC, I

      WRITE(JICHAR(1),'(F10.5)')JLOC
      WRITE(JICHAR(2),'(F10.5)')ILOC
      print*, 'JLOC = ',JLOC,' ILOC = ',ILOC,' DPREC = ',DPREC
      CLEN_OLD=1
      DO I=1,2
         NSPACE=1
 5       IF(JICHAR(I)(NSPACE:NSPACE).EQ.' ')THEN
            NSPACE=NSPACE+1
            GOTO 5
         ENDIF
         CLEN=CLEN_OLD+11-NSPACE-5+DPREC
         EXTEN(CLEN_OLD:CLEN)=JICHAR(I)(NSPACE:5+DPREC)
      print*,i,clen,clen_old,nspace,dprec,exten
         IF(I.EQ.1)THEN
            EXTEN(CLEN:CLEN)='_'
         ENDIF
         CLEN_OLD=CLEN+1
      END DO
      pause
      CLEN=CLEN-1
      
      RETURN
      END
