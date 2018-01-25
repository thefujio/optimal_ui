SUBROUTINE movavg(invec,outvec,pd)
!Simple moving average of 1-D real*8 vector with odd period = pd

implicit none

real*8, dimension(:), intent(in) :: invec
integer, intent(in) ::pd
real*8, dimension(:), intent(out) :: outvec
real*8 :: sum
integer :: i,j,jj,len,hpd

len = size(invec)

hpd = floor(real(pd,8)/2.0d0)

do i=1,len
  sum = 0.0d0
  do j = 1,hpd
    sum = sum + invec(min(len,i+j))
  end do
  do jj = 1,hpd
    sum = sum + invec(max(1,i-jj))
  end do
    sum = sum + invec(i)
    outvec(i) = sum / real(pd,8)
end do

END SUBROUTINE movavg
