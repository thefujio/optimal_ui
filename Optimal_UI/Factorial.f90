RECURSIVE SUBROUTINE Factorial(N,Result)
  INTEGER, INTENT(IN):: N
  INTEGER, INTENT(INOUT):: Result
  
  IF (N>0) THEN
    CALL Factorial(N-1,Result)
    Result = Result*N
  ELSE
    Result = 1
  END IF
END SUBROUTINE Factorial
