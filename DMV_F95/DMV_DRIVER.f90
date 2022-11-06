!******************************************************************************
!                                                                              
!  PROGRAM: DMV_DRIVER                                                         
!                                                                              
!  PURPOSE:  DRIVER FOR DMV10 FOR THE MOTION OF A RIGID BODY                   
!            DMV10 : PREPROCESSED DISCRETE MOSER-VESELOV ALGORITHMS OF ORDER=10
!                                                                              
!******************************************************************************

    program DMV_DRIVER
    use mod_dmv
    implicit none

    ! Variables
    type(mmoi) :: PAR
    REAL(rk) :: Q(4), AM(3)
    REAL(rk) :: H, XEND, HAM0, HAM1
    INTEGER :: NSTEP
    REAL :: TIME
    INTEGER(int64) :: tic,toc,rate
    
! --- INITIAL VALUES
    
! ANGULAR MOMENTUM
        AM = [ 1.8_rk, 0.4_rk, -0.9_rk]
        
! QUATERNION
        Q = [1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
        
! --- MOMENTS OF INERTIA
        PAR = MMOI(0.6_rk,0.8_rk,1.0_rk)

        H = 1e-6_rk
        XEND = 10_rk
        NSTEP=XEND/H
        H=XEND/NSTEP
        WRITE (*,*) 'XEND=',XEND,' H=',H,' NSTEP=',NSTEP

        WRITE (*,*) '--- INITIAL CONDITION'
        WRITE (*,*) '  AM    ', REAL(AM)
        WRITE (*,*) '  QQ    ', REAL(Q)
        HAM0 = HAMIL(Q,AM,PAR)
        WRITE (*,*) '  HAM   ', REAL(HAM0)
        CALL SYSTEM_CLOCK(tic,rate)

        CALL DMV10(AM,Q,POTENP,H,NSTEP,PAR)

        CALL SYSTEM_CLOCK(toc,rate)
        TIME = (toc-tic)/float(rate)
        HAM1 = HAMIL(Q,AM,PAR)
        WRITE (*,*) '--- SOLUTION AT ENDPOINT'
        WRITE (*,*) '  AM    ', REAL(AM)
        WRITE (*,*) '  QQ    ', REAL(Q)        
        WRITE (*,*) '  HAM   ', REAL(HAM1)
        WRITE (*,*) ''
        WRITE (*,*) 'ERR HAM=',REAL(HAM1/HAM0-1.0_rk),' TIME=',time &
                   ,' MFPS=',REAL(NSTEP/time/1e6)

        STOP
    contains
    
        PURE FUNCTION POTENP(Q,PAR) RESULT(POTP)
        import :: rk
        REAL(rk), INTENT(IN) :: Q(4)
        type(mmoi), intent(in) :: PAR
        REAL(rk) :: POTP(3)
        !POTP(1)=-2*(Q(1)*Q(2)+Q(3)*Q(4))
        !POTP(2)=-2*(Q(1)*Q(3)-Q(2)*Q(4))
        POTP = [0.0_rk, 0.0_rk, 0.0_rk]
        RETURN
        END FUNCTION
        
        PURE FUNCTION HAMIL(Q,AM,PAR) RESULT(HAM)
        import :: rk
        REAL(rk), INTENT(IN) :: Q(4), AM(3)
        type(mmoi), intent(in) :: PAR
        REAL(rk) :: HAM
        REAL(rk) :: POT
        HAM=AM(1)**2/PAR%AI1+AM(2)**2/PAR%AI2+AM(3)**2/PAR%AI3
        !POT=Q(1)**2-Q(2)**2-Q(3)**2+Q(4)**2
        POT=0.0_rk    
        HAM=HAM/2+POT
        RETURN
        END FUNCTION
        
    
    end program

