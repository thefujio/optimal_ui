MODULE UTILITY
  !*************************************************************************
  !  UTILITY.f90
  !  
  !   
  !  PURPOSE:
  !   Makes available all the functions related to the utility function and
  !   its derivatives
  !
  !
  !  Program written by: M. Gervais
    !   Date:    Jan 2016
    !   Updated: Feb 2016
  !*************************************************************************
  implicit none
  private
  public:: Ufunc,Um1

CONTAINS
  
  DOUBLE PRECISION FUNCTION Ufunc(c)
    !****************************************************************************
    !  Ufunc.f90
    !  
    !  PURPOSE:
    !   Computes utility for an individual consuming c
    !
    !
    !  Program written by: M. Gervais	
    !   Date:    Jan 2016
    !   Updated: Feb 2016
    !****************************************************************************
    USE PARAM
    implicit none
    
    !dummy arguments declarations
    double precision, intent(in):: c
    
    !Local variables declarations
    
    if (c.LE.zero) then
      Ufunc = hell
      RETURN
    end if
    
    if (dabs(sigma-one)<tol) then
      Ufunc = dlog(c)
    else
      Ufunc = (c**(one-sigma)-one)/(one-sigma)
    end if
    
    RETURN
  END FUNCTION Ufunc
  
  
  DOUBLE PRECISION FUNCTION Um1(utility)
    !****************************************************************************
    !  Uc.f90
    !  
    !  PURPOSE:
    !   Computes the implied consumption if utility is equal to u.
    !
    !
    !  Program written by: M. Gervais
    !   Date:    Feb 2016
    !   Updated: Feb 2016
    !****************************************************************************
    USE PARAM
    implicit none
    
    !dummy arguments declarations
    double precision, intent(in):: utility
    
    if (dabs(sigma-one)<tol) then
      Um1 = exp(utility)
    else
      Um1 = (utility*(one-sigma)+one)**(one/(one-sigma))
    end if
    
    if (Um1.LE.zero) then
      Um1 = hell
      RETURN
    end if
    
    RETURN
  END FUNCTION Um1
END MODULE UTILITY
