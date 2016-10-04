program combinations
 
  implicit none
  integer, parameter :: m_max = 3
  integer, parameter :: n_max = 10
  integer, dimension (m_max) :: comb
  character (*), parameter :: fmt = '(i0' // repeat (', 1x, i0', m_max - 1) // ')'
  integer:: n,i,j
  integer,allocatable:: P(:,:)
  
  n=3
  allocate( P(product((/(i, i=1,n)/)),  n) ) 
  call permutate( (/(i, i=1,n)/),  P ) 
  do j=1,size(P,1) 
    print *, P(j,:) 
  end do 

! 
!  call gen (1)
 
contains
 
  recursive subroutine gen (m)
 
    implicit none
    integer, intent (in) :: m
    integer :: n
 
    if (m > m_max) then
      write (*, fmt) comb
    else
      do n = 1, n_max
        if ((m == 1) .or. (n > comb (m - 1))) then
          comb (m) = n
          call gen (m + 1)
        end if
      end do
    end if
 
  end subroutine gen
  
  recursive subroutine permutate(E, P) 
    integer, intent(in)  :: E(:)       ! array of objects 
    integer, intent(out) :: P(:,:)     ! permutations of E 
    integer  :: N, Nfac, i, k, S(size(P,1)/size(E), size(E)-1) 
    N = size(E); Nfac = size(P,1); 
    do i=1,N                           ! cases with E(i) in front 
      if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S) 
      forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/) 
    end do 
  end subroutine permutate 
  
 
end program combinations
