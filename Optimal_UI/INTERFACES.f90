MODULE INTERFACES
  !****************************************************************************
  !  INTERFACE.f90
  !  
  !  PURPOSE:
  !   Sets all interfaces
  !  
  !  
  !  Program written by: J. Fisher and M. Gervais
  !   Date:    April 2005
  !   Updated: June 2008
  !****************************************************************************
  implicit none
  
  INTERFACE
    DOUBLE PRECISION FUNCTION c_l_foc(l)
      USE PARAM
      USE UTILITY
      USE PRODUCTION
      implicit none
      double precision:: l
    END FUNCTION c_l_foc
    
    DOUBLE PRECISION FUNCTION c0_l0_foc(l)
      USE PARAM
      USE UTILITY
      USE PRODUCTION
      implicit none
      double precision:: l
    END FUNCTION c0_l0_foc
    
    DOUBLE PRECISION FUNCTION css_lss_foc(l)
      USE PARAM
      USE UTILITY
      USE PRODUCTION
      implicit none
      double precision:: l
    END FUNCTION css_lss_foc
    
    DOUBLE PRECISION FUNCTION euler_eq(kp)
      USE PARAM
      USE UTILITY
      USE PRODUCTION
      implicit none
      double precision:: kp
    END FUNCTION euler_eq
  END INTERFACE
END MODULE INTERFACES