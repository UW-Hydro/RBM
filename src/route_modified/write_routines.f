c SUBROUTINES RELATED TO WRITING
c write_data()
c write_month()
c

      SUBROUTINE WRITE_DATA
c
c  Add Qin - JRY - 9/30/2009
c 
     & (write_flow,rbm10,FLOW, Qin, DAYS,GRID_CELL
     & ,FLOWOUT,HEATPATH,clen
     & ,Flow_Cells,Force_Cells,Nflow,Nflow_pd,Nheat,Nheat_pd 
     & ,nday,IDAY,IMONTH,IYEAR,skip_MET
     & ,aa_d,bb_d,aa_w,bb_w)
c 
      IMPLICIT NONE

      INTEGER DAYS
c
c  Add Qin - JRY - 9/30/2009
c
      REAL    FLOW(DAYS),Qin(days),Dmmy(7),Heat_data(7)
      INTEGER IDAY(DAYS), IMONTH(DAYS), IYEAR(DAYS)
      INTEGER I, ii, CLEN, FLEN, Flow_Cells, Force_Cells
     &       ,nday,Nflow,Nflow_pd, navg, Nheat, Nheat_pd, nrec
     &       ,skip_MET,nyy,nmm,ndd,nhh
     &       ,Ntest
      logical rbm10,write_flow,write_header
      CHARACTER*200 GRID_CELL 
      REAL    FACTOR_SUM
     &    ,flowavg,flowin,depth,width,vel,heat_pd
     &    ,yy,mm,dd,hh
      real    aa_d,bb_d,aa_w,bb_w
c 
      character*200 header
      CHARACTER*200 OUTPATH
      character*200 flowout,flowpath,heatout,heatpath
      if (rbm10) then
         heat_pd=Nheat_pd
c
c     Open forcing file from VIC
c 
         write(*,'(A)') 'heat path ',TRIM(heatpath)//TRIM(GRID_CELL)
         open(25,FILE=TRIM(heatpath)//TRIM(GRID_CELL),status='old')
c
c     Read header on forcing file
c 
        do ii=1,6
          read(25,'(a)') header
        end do
c
c     Skip records if necessary
        if (skip_MET.gt.0) then
           do ii=1,skip_MET
              read(25,*)
           end do
         end if  
         do i=1,nday
           if (write_flow) then
              nrec=Flow_Cells*(i-1)+Nflow
              flowin=flow(i)-Qin(i)
              if(flowin.lt.5.0) flowin=5.0
              flowavg=0.5*(flowin+flow(i))
              depth=aa_d*(flowavg**bb_d)
              width=aa_w*(flowavg**bb_w)
              vel=flowavg/(depth*width)
c      
c     Write the flow for this grid cell to the RBM direct access file
c 
             write(15,'(2i5,2f10.1,2f6.1,f7.1,f6.2)',rec=nrec)
     &              i,Nheat,flowin,flow(i),Qin(i),depth,width,vel
           end if
             do navg=1,7
                Heat_data(navg)=0.0
             end do
             do navg=1,Nheat_pd
c
c            Read full_data output from VIC runs and average
c            year,month,day,air_temp,vp,short_wave,long_wave,
c            rho,press,wind
c 
               read(25,*) dmmy
               Heat_data(1)=Heat_data(1)+dmmy(1)/heat_pd
               Heat_data(2)=Heat_data(2)+10.*dmmy(2)/heat_pd
               Heat_data(3)=Heat_data(3)+2.388e-4*dmmy(3)/heat_pd
               Heat_data(4)=Heat_data(4)+2.388e-4*dmmy(4)/heat_pd
               Heat_data(5)=Heat_data(5)+dmmy(5)/heat_pd
               Heat_data(6)=Heat_data(6)+10.*dmmy(6)/heat_pd
               Heat_data(7)=Heat_data(7)+dmmy(7)/heat_pd
            end do
            nrec=Force_Cells*(i-1)+Nheat
c
c         Write the heat budget data for this grid cell
c         to the RB direct access file
c         
            write(16,'(i5,2f6.1,2f7.4,f6.3,f7.1,f5.1)',rec=nrec)
     &           Nheat,Heat_data
         end do
         close(25)
      else 
         OPEN(30, FILE = TRIM(FLOWOUT)//TRIM(GRID_CELL)//'.day'
     &          ,status='unknown')
         DO I = 1,DAYS
c
c  Add Qin - JRY- 9/30/2009
c 
            WRITE(30,*) IYEAR(I),IMONTH(I),IDAY(I),FLOW(I),Qin(i)
         END DO
      end if
      CLOSE(30)
      RETURN
      END
