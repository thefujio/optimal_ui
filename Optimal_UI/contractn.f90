SUBROUTINE contract(cont,n,k,nn)
  !****************************************************************************
  !  contractn.f90
  !  Computes all possible contracts (indexes) given k states and n submarkets
  !  Larry Warren
  !  September 2016
  !****************************************************************************
  
  implicit none
  !dummy argument declarations
  integer, intent(in):: n,k,nn
  integer, intent(out)::cont(nn,k)
  
  !Local variable declarations
  integer:: row,col,ii
  integer:: indices(k)
  integer:: P(nn,k)
  indices(:) = 1
    do row = 1,nn
	  do col = 1,k
	    P(row, col) = indices(col)
      end do
	  indices(k) = indices(k) + 1
	  ii = k
	  do while (ii>1 .and. indices(ii) > n)
	    indices(ii) = 1
	    indices(ii-1) = indices(ii -1) + 1
	    ii = ii-1
	  end do
    enddo
	
  cont = P
  
  RETURN
END SUBROUTINE contract