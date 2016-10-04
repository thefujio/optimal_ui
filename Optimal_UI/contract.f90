SUBROUTINE contract(cont,n,k,nn)
  !****************************************************************************
  !  contract.f90
  !  
  !  PURPOSE:
  !   Computes all possible contracts given k (number of states) and n (number 
  !   of submarkets)
  !
  !
  !  Program written by: M. Gervais	
  !   Date:    Feb 2016
  !   Updated: Feb 2016
  !****************************************************************************
  implicit none
  
  !dummy arguments declarations
  integer, intent(in):: n,k,nn
  integer, intent(out):: cont(nn,k)
  
  !Local variables declarations
  integer:: i,j,i1,i2,i3,i4
  integer:: P(nn,k)
  
  
  j=1
  if (k==1) then
    P(:,1) = (/(i, i=1,n)/)
  else if (k==2) then
    do i1=1,n
      do i2=1,n
        P(j,1) = i1
        P(j,2) = i2
        j=j+1
      end do
    end do
  else if (k==3) then
    do i1=1,n
      do i2=1,n
        do i3=1,n
          P(j,1) = i1
          P(j,2) = i2
          P(j,3) = i3
          j=j+1
        end do
      end do
    end do
  else if (k==4) then
    do i1=1,n
      do i2=1,n
        do i3=1,n
          do i4=1,n
            P(j,1) = i1
            P(j,2) = i2
            P(j,3) = i3
            P(j,4) = i4
            j=j+1
          end do
        end do
      end do
    end do
  end if
  
  cont = P
  
  RETURN
END SUBROUTINE contract
