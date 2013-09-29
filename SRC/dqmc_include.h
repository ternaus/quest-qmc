!
! DEFINITION OF SOME PRE_PROCESSOR

!
! ==================================================  
!   CKB means checker-board method for sparse matrix B.
!   If it is defined, the program will use check-board 
!   method in B matrix multiplication.
!   Otherwise, a dense matrix B will be formed 
!   explicitly.
! ==================================================  
!
#ifdef  _CKB
#define _DQMC_MATB DQMC_CKB
#else
#define _DQMC_MATB DQMC_MATB
#endif

