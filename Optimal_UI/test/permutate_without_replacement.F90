program combinations
 
  implicit none
  integer, parameter :: k = 2
  integer, parameter :: n = 4
  integer::i
  integer:: kfac,nfac
  integer, dimension (k) :: comb
  character (*), parameter :: fmt = '(i0' // repeat (', 1x, i0', k - 1) // ')'
  integer,allocatable:: P(:,:)
  
  kfac = product((/(i, i=1,k)/))
  nfac = product((/(i, i=1,n)/))
  allocate( P(nfac/(kfac*(nfac-kfac)),  k) )
 
  call gen (1,n,k,P)
  
  do i=1,6
    write (*, fmt) P(i,:)
  end do
 
contains
 
  recursive subroutine gen (m,n,k,P)
 
    implicit none
    integer, intent (in) :: m,n,k
    integer, intent (out):: P(:,:)
    integer:: i,row
    
    row = 1
    if (m > k) then
      write (*, fmt) comb
      P(row,:) = comb
      write (*, fmt) P(i,:)
      row = row + 1
    else
      do i = 1, n
        if ((m == 1) .or. (i > comb (m - 1))) then
          comb (m) = i
          call gen (m+1,n,k,P)
        end if
      end do
    end if
 
  end subroutine gen
 
end program combinations
