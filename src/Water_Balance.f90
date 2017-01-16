SUBROUTINE Water_Balance
!
USE Block_Hydro
USE Block_Network
!
IMPLICIT NONE
!
integer:: nc,ncnf,nr,nrc,ntrib,nntrib,tot_trib_in,trib_seg_in
real   :: Q_start,Q_sum_trib
!
do nr = 1,nreach
  nrc=reach_cell(nr,1)
  Q_in(nrc) = Q_out(nrc)
  Q_diff(nrc)=0.0
  nc=1
  write(27,*) 'Water Balance - ',nr,nc,nrc,Q_diff(nrc),Q_out(nrc),Q_in(nrc)  &
                                  ,Q_sum_trib
  do nc = 2,no_cells(nr)-1
    nrc = reach_cell(nr,nc)
    Q_in(nrc) = Q_out(nrc-1)
    ntrib=no_tribs(nrc)
    Q_sum_trib = 0.0
    if (ntrib .gt. 0) then
      do nntrib  = 1,ntrib
        ncnf = conflnce(nrc,nntrib)
        Q_sum_trib = Q_sum_trib + Q_out(ncnf)
      end do
    end if
    Q_diff(nrc) = Q_out(nrc) - Q_in(nrc) - Q_sum_trib
    write(27,*) 'Water Balance - ',nr,nc,nrc,Q_diff(nrc),Q_out(nrc),Q_in(nrc)  &
                                  ,Q_sum_trib
  end do
end do
end SUBROUTINE Water_Balance
