    module mod_dmv
    !-----------------------------------------------------------------------
    !            PREPROCESSED DISCRETE MOSER-VESELOV ALGORITHMS             
    !-----------------------------------------------------------------------
    !  PREPROCESSED DISCRETE MOSER-VESELOV ALGORITHM OF ORDERS 2-10 FOR THE 
    !  NUMERICAL SOLUTION OF THE EQUATIONS OF MOTION OF THE FREE RIGID BODY. 
    !  ORTHOGONAL MATRICES ARE REPRESENTED BY QUATERNIONS. 
    !  THE CODE IS READY TO INCLUDE AN EXERNAL POTENTIAL (SYMMETRIC STRANG  
    !  SPLITTING OF ORDER 2). 
    !-----------------------------------------------------------------------
    !     AUTHORS: ERNST HAIRER (1) AND GILLES VILMART (1)(2)
    !                (1) UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
    !                    2-4 RUE DU LIEVRE, CASE POSTALE 64
    !                    CH-1211 GENEVE 4, SWITZERLAND 
    !                (2) IRISA/INRIA RENNES, PROJET IPSO
    !                    CAMPUS DE BEAULIEU, F-35042 RENNES CEDEX, FRANCE
    !     E-MAILS: Ernst.Hairer@math.unige.ch
    !              Gilles.Vilmart@math.unige.ch
    !
    !  THIS CODE IS DESCRIBED IN:
    !    E. HAIRER AND G. VILMART, PREPROCESSED DISCRETE MOSER-VESELOV
    !    ALGORITHM FOR THE FULL DYNAMICS OF A RIGID BODY
    !    J. Phys. A: Math. Gen. 39 (2006) 13225-13235. 
    !    http://stacks.iop.org/0305-4470/39/13225
    !    [pdf] => https://www.unige.ch/~hairer/preprints/predmv.pdf
    !
    !     VERSION: AUGUST 29, 2006
    !     (latest correction of a small bug: December 28, 2007)
    !-----------------------------------------------------------------------
        
    
    ! Standardized Fortran 2008 kind definitions
    use, intrinsic :: ISO_FORTRAN_ENV, only : &
        INT8, INT16, INT32, INT64, &       
        REAL32, REAL64, REAL128
    !IMPLICIT REAL(real64) (A-H,O-Z)
    IMPLICIT NONE    
    integer, parameter :: rk = REAL64
    
    public
    
    type mmoi
      ! THREE MOMENTS OF INERTIA OF THE RIGID BODY
      real(rk) :: AI1, AI2, AI3
    end type
    
    INTERFACE
    PURE FUNCTION POTENTIAL(Q,PAR) RESULT(POTP)
    IMPORT
    REAL(RK), INTENT(IN) :: Q(4)
    type(mmoi), intent(in) :: PAR
    REAL(RK) :: POTP(3)
    END FUNCTION
    END INTERFACE

    real(rk), parameter :: tiny = 1E-15_rk
    contains
    
    ! SUBROUTINE DMVn(AM,Q,POTENP,H,NSTEP,PAR)
    !
    !     INPUT PARAMETERS
    !     ----------------
    !     AM(I)      INITIAL ANGULAR MOMENTUM (I=1,2,3)
    !     Q(I)       INITIAL QUATERNION FOR ORTHOGONAL MATRIX (I=1,2,3,4) 
    !     H          STEP SIZE
    !     NSTEP      NUMBER OF STEPS
    !     POTENP     NAME (EXTERNAL) OF SUBROUTINE FOR AN EXTERNAL POTENTIAL
    !                   FUNCTION POTENP(Q,PAR) RESULT(POTP)
    !                   DIMENSION Q(4),POTP(3)
    !                   POTP(1)=...   ETC.
    !     PAR        MASS MOMENT OF INERTIA PARAMETERS, user-type "mmoi"
    !
    !     OUTPUT PARAMETERS
    !     -----------------
    !     AM(I)       SOLUTION (ANGULAR MOMENTUM) AT ENDPOINT
    !     Q(I)        SOLUTION (QUATERNION) AT ENDPOINT
        
    PURE SUBROUTINE DMV2(AM,Q,POTENP,H,NSTEP,PAR)
    !--- DISCRETE MOSER-VESELOV ALGORITHM OF ORDER 2
    !--- ORTHOGONAL MATRICES ARE REPRESENTED BY QUATERNIONS
    
      ! Arguments
      procedure(POTENTIAL) :: POTENP
      REAL(rk), intent(inout) :: Q(4),AM(3)
      REAL(rk), intent(in) :: H
      INTEGER, intent(in) :: NSTEP
      type(mmoi), intent(in) :: PAR
    
      ! Subroutine Constants
      integer, PARAMETER :: MSPLIT=1
      
      ! Local Variables
      REAL(rk) :: POTP(3)                                                       ! POTENTIAL ENERGY
      REAL(rk) :: EPS, HA, HC, ERR                                      ! TIME STEP RELATED
      REAL(rk) :: AI1, AI2, AI3, AM1, AM2, AM3                                  ! MMOI VALUES
      REAL(rk) :: FAD1, FAD2, FAD3, FAC1, FAC2, FAC3
      REAL(rk) :: AM1I, AM2I, AM3I, CM1, CM2, CM3, CM1B, CM2B, CM3B, CALPHA
      REAL(rk) :: Q0, Q1, Q2, Q3, QUAT
      
      INTEGER :: ISTEP, I
      
      EPS=ABS(tiny*H)
      HA=H/2
      AI1=PAR%AI1
      AI2=PAR%AI2
      AI3=PAR%AI3 
      HC=4/H
      FAD1=HC*(AI2-AI3)
      FAD2=HC*(AI3-AI1)
      FAD3=HC*(AI1-AI2)
      FAC1=(AI2-AI3)/AI1
      FAC2=(AI3-AI1)/AI2
      FAC3=(AI1-AI2)/AI3
!---
      POTP = POTENP(Q,PAR)
      AM1=AM(1)-HA*POTP(1)
      AM2=AM(2)-HA*POTP(2)
      AM3=AM(3)-HA*POTP(3)

      DO ISTEP=1,NSTEP
! ---      SOLVE FOR INTERNAL STAGE
        AM1I=AM1*HA/AI1
        AM2I=AM2*HA/AI2
        AM3I=AM3*HA/AI3
        CM1=AM1I+FAC1*AM2I*AM3I
        CM2=AM2I+FAC2*CM1*AM3I
        CM3=AM3I+FAC3*CM1*CM2
        DO I=1,50
          CM1B=CM1
          CM2B=CM2
          CM3B=CM3
          CALPHA=1+CM1**2+CM2**2+CM3**2
          CM1=CALPHA*AM1I+FAC1*CM2*CM3
          CM2=CALPHA*AM2I+FAC2*CM1*CM3
          CM3=CALPHA*AM3I+FAC3*CM1*CM2
          ERR=ABS(CM1B-CM1)+ABS(CM2B-CM2)+ABS(CM3B-CM3)
          IF( ERR < EPS ) EXIT
        END DO
! ---      UPDATE Q
        Q0=Q(1)
        Q1=Q(2)
        Q2=Q(3)
        Q3=Q(4)
        Q(1)=Q0-CM1*Q1-CM2*Q2-CM3*Q3
        Q(2)=Q1+CM1*Q0+CM3*Q2-CM2*Q3
        Q(3)=Q2+CM2*Q0+CM1*Q3-CM3*Q1
        Q(4)=Q3+CM3*Q0+CM2*Q1-CM1*Q2
! ---       PROJECTION
        QUAT=1/SQRT(Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2)
        Q(1)=Q(1)*QUAT
        Q(2)=Q(2)*QUAT
        Q(3)=Q(3)*QUAT
        Q(4)=Q(4)*QUAT
! ---       UPDATE M
        POTP = POTENP(Q,PAR)
        AM1=AM1-H*POTP(1)+FAD1*CM2*CM3/CALPHA
        AM2=AM2-H*POTP(2)+FAD2*CM1*CM3/CALPHA
        AM3=AM3-H*POTP(3)+FAD3*CM1*CM2/CALPHA
      END DO
      AM(1)=AM1+HA*POTP(1)
      AM(2)=AM2+HA*POTP(2)
      AM(3)=AM3+HA*POTP(3)
      RETURN
      
    END SUBROUTINE
    
    PURE SUBROUTINE DMV4(AM,Q,POTENP,H,NSTEP,PAR)
    !--- DISCRETE MOSER-VESELOV ALGORITHM OF ORDER 4
    !--- ORTHOGONAL MATRICES ARE REPRESENTED BY QUATERNIONS
    
      ! Arguments
      procedure(POTENTIAL) :: POTENP
      REAL(rk), intent(inout) :: Q(4),AM(3)
      REAL(rk), intent(in) :: H
      INTEGER, intent(in) :: NSTEP
      type(mmoi), intent(in) :: PAR
    
      ! Subroutine Constants
      integer, PARAMETER :: MSPLIT=1
      
      ! Local Variables
      REAL(rk) :: POTP(3)                                                       ! POTENTIAL ENERGY
      REAL(rk) :: EPS, HA, HC, H2, ERR                                      ! TIME STEP RELATED
      REAL(rk) :: AI1, AI2, AI3, AM1, AM2, AM3                                  ! MMOI VALUES
      REAL(rk) :: FAD1, FAD2, FAD3, FAC1, FAC2, FAC3
      REAL(rk) :: AM1I, AM2I, AM3I, CM1, CM2, CM3, CM1B, CM2B, CM3B, CALPHA
      REAL(rk) :: Q0, Q1, Q2, Q3, QUAT
      REAL(rk) :: XI1, XDET, XS1
      REAL(rk) :: D3C1, D3C2
      REAL(rk) :: HAM0, S3C1, S3C2
      REAL(rk) :: ANOR0
      REAL(rk) :: CSS, CDD, AI1MODI, AI2MODI, AI3MODI, AI1MOD, AI2MOD, AI3MOD
      
      INTEGER :: ISTEP, I
      
      EPS=ABS(tiny*H)
      HA=H/2
      H2 = H**2
      AI1=PAR%AI1
      AI2=PAR%AI2
      AI3=PAR%AI3 
      HC=4/H

      XI1=1/AI1+1/AI2+1/AI3
      XDET=AI1*AI2*AI3
      XS1=AI1+AI2+AI3
      S3C1=XS1/(6*XDET)
      S3C2=-XI1/3
      D3C1=-1/(3*XDET)
      D3C2=S3C1

      POTP = POTENP(Q,PAR)
      AM1=AM(1)-HA*POTP(1)
      AM2=AM(2)-HA*POTP(2)
      AM3=AM(3)-HA*POTP(3)

      DO ISTEP=1,NSTEP
        HAM0=0.5_rk*(AM1**2/AI1+AM2**2/AI2+AM3**2/AI3)
        ANOR0=0.5_rk*(AM1**2+AM2**2+AM3**2)
        CSS=1+H2*(S3C1*ANOR0+S3C2*HAM0)
        CDD=H2*(D3C1*ANOR0+D3C2*HAM0)
        AI1MODI=CSS/AI1+CDD
        AI2MODI=CSS/AI2+CDD
        AI3MODI=CSS/AI3+CDD
        AI1MOD=1/AI1MODI
        AI2MOD=1/AI2MODI
        AI3MOD=1/AI3MODI
        FAD1=AI2MOD-AI3MOD
        FAD2=AI3MOD-AI1MOD
        FAD3=AI1MOD-AI2MOD
        FAC1=FAD1*AI1MODI
        FAC2=FAD2*AI2MODI
        FAC3=FAD3*AI3MODI
! ---      SOLVE FOR INTERNAL STAGE
        AM1I=AM1*HA*AI1MODI
        AM2I=AM2*HA*AI2MODI
        AM3I=AM3*HA*AI3MODI
        CM1=AM1I+FAC1*AM2I*AM3I
        CM2=AM2I+FAC2*CM1*AM3I
        CM3=AM3I+FAC3*CM1*CM2
        DO I=1,50
          CM1B=CM1
          CM2B=CM2
          CM3B=CM3
          CALPHA=1+CM1**2+CM2**2+CM3**2
          CM1=CALPHA*AM1I+FAC1*CM2*CM3
          CM2=CALPHA*AM2I+FAC2*CM1*CM3
          CM3=CALPHA*AM3I+FAC3*CM1*CM2
          ERR=ABS(CM1B-CM1)+ABS(CM2B-CM2)+ABS(CM3B-CM3)
          IF( ERR < EPS ) EXIT
        END DO
! ---      UPDATE Q
        Q0=Q(1)
        Q1=Q(2)
        Q2=Q(3)
        Q3=Q(4)
        Q(1)=Q0-CM1*Q1-CM2*Q2-CM3*Q3
        Q(2)=Q1+CM1*Q0+CM3*Q2-CM2*Q3
        Q(3)=Q2+CM2*Q0+CM1*Q3-CM3*Q1
        Q(4)=Q3+CM3*Q0+CM2*Q1-CM1*Q2
! ---       PROJECTION
        QUAT=1/SQRT(Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2)
        Q(1)=Q(1)*QUAT
        Q(2)=Q(2)*QUAT
        Q(3)=Q(3)*QUAT
        Q(4)=Q(4)*QUAT
! ---       UPDATE M
        POTP = POTENP(Q,PAR)
        CALPHA=HC/CALPHA
        AM1=AM1-H*POTP(1)+FAD1*CM2*CM3*CALPHA
        AM2=AM2-H*POTP(2)+FAD2*CM1*CM3*CALPHA
        AM3=AM3-H*POTP(3)+FAD3*CM1*CM2*CALPHA
      END DO
      AM(1)=AM1+HA*POTP(1)
      AM(2)=AM2+HA*POTP(2)
      AM(3)=AM3+HA*POTP(3)
      RETURN      
    END SUBROUTINE
    
    PURE SUBROUTINE DMV6(AM,Q,POTENP,H,NSTEP,PAR)
    !--- DISCRETE MOSER-VESELOV ALGORITHM OF ORDER 6
    !--- ORTHOGONAL MATRICES ARE REPRESENTED BY QUATERNIONS
    
      ! Arguments
      procedure(POTENTIAL) :: POTENP
      REAL(rk), intent(inout) :: Q(4),AM(3)
      REAL(rk), intent(in) :: H
      INTEGER, intent(in) :: NSTEP
      type(mmoi), intent(in) :: PAR
    
      ! Subroutine Constants
      integer, PARAMETER :: MSPLIT=1
      
      ! Local Variables
      REAL(rk) :: POTP(3)                                                       ! POTENTIAL ENERGY
      REAL(rk) :: EPS, HA, HC, H2, ERR                                      ! TIME STEP RELATED
      REAL(rk) :: AI1, AI2, AI3, AM1, AM2, AM3                                  ! MMOI VALUES
      REAL(rk) :: FAD1, FAD2, FAD3, FAC1, FAC2, FAC3
      REAL(rk) :: AM1I, AM2I, AM3I, CM1, CM2, CM3, CM1B, CM2B, CM3B, CALPHA
      REAL(rk) :: Q0, Q1, Q2, Q3, QUAT, XI1, XI2, XDET, XS1, XS2, XSI, XSJ
      REAL(rk) :: S3C1, S3C2, D3C1, D3C2, S5C1, S5C2, S5C3, D5C1, D5C2, D5C3
      REAL(rk) :: HAM0, ANOR0, HAM2, ANOR2, ANORHAM, CSS, CDD
      REAL(rk) :: AI1MODI, AI2MODI, AI3MODI, AI1MOD, AI2MOD, AI3MOD
      
      INTEGER :: ISTEP, I
      
      EPS=ABS(tiny*H)
      
      HA=H/2
      HC=4/H
      H2 = H**2
      
      AI1=PAR%AI1
      AI2=PAR%AI2
      AI3=PAR%AI3 

      XI1=1/AI1+1/AI2+1/AI3
      XI2=1/AI1**2+1/AI2**2+1/AI3**2
      XDET=AI1*AI2*AI3
      XS1=AI1+AI2+AI3
      XS2=AI1**2+AI2**2+AI3**2
      XSI=(AI1+AI2)/AI3+(AI2+AI3)/AI1+(AI3+AI1)/AI2
      XSJ=AI1/AI2/AI3+AI2/AI3/AI1+AI3/AI1/AI2

      S3C1=XS1/(6*XDET)
      S3C2=-XI1/3
      D3C1=-1/(3*XDET)
      D3C2=S3C1

      S5C1=(XS2/XDET-XI1)/(30*XDET)
      S5C2=(1-XSI)/(30*XDET)
      S5C3=(3*XS1/XDET+2*XI2)/60
      D5C1=-XS1/XDET**2/60
      D5C2=XI1/XDET/10-XS2/XDET**2/60
      D5C3=-(9+XSI)/XDET/60

      POTP = POTENP(Q,PAR)
      AM1=AM(1)-HA*POTP(1)
      AM2=AM(2)-HA*POTP(2)
      AM3=AM(3)-HA*POTP(3)

      DO ISTEP=1,NSTEP
        HAM0=0.5_rk*H2*(AM1**2/AI1+AM2**2/AI2+AM3**2/AI3)
        ANOR0=0.5_rk*H2*(AM1**2+AM2**2+AM3**2)
        HAM2=HAM0**2
        ANOR2=ANOR0**2
        ANORHAM=HAM0*ANOR0
        CSS=1+(S3C1*ANOR0+S3C2*HAM0) &
           +(S5C1*ANOR2+S5C2*ANORHAM+S5C3*HAM2)
        CDD=(D3C1*ANOR0+D3C2*HAM0) &
           +(D5C1*ANOR2+D5C2*ANORHAM+D5C3*HAM2)
        AI1MODI=CSS/AI1+CDD
        AI2MODI=CSS/AI2+CDD
        AI3MODI=CSS/AI3+CDD
        AI1MOD=1/AI1MODI
        AI2MOD=1/AI2MODI
        AI3MOD=1/AI3MODI
        FAD1=AI2MOD-AI3MOD
        FAD2=AI3MOD-AI1MOD
        FAD3=AI1MOD-AI2MOD
        FAC1=FAD1*AI1MODI
        FAC2=FAD2*AI2MODI
        FAC3=FAD3*AI3MODI
! ---      SOLVE FOR INTERNAL STAGE
        AM1I=AM1*HA*AI1MODI
        AM2I=AM2*HA*AI2MODI
        AM3I=AM3*HA*AI3MODI
        CM1=AM1I+FAC1*AM2I*AM3I
        CM2=AM2I+FAC2*CM1*AM3I
        CM3=AM3I+FAC3*CM1*CM2
        DO I=1,50
          CM1B=CM1
          CM2B=CM2
          CM3B=CM3
          CALPHA=1+CM1**2+CM2**2+CM3**2
          CM1=CALPHA*AM1I+FAC1*CM2*CM3
          CM2=CALPHA*AM2I+FAC2*CM1*CM3
          CM3=CALPHA*AM3I+FAC3*CM1*CM2
          ERR=ABS(CM1B-CM1)+ABS(CM2B-CM2)+ABS(CM3B-CM3)
          IF( ERR < EPS ) EXIT
        END DO
! ---      UPDATE Q
        Q0=Q(1)
        Q1=Q(2)
        Q2=Q(3)
        Q3=Q(4)
        Q(1)=Q0-CM1*Q1-CM2*Q2-CM3*Q3
        Q(2)=Q1+CM1*Q0+CM3*Q2-CM2*Q3
        Q(3)=Q2+CM2*Q0+CM1*Q3-CM3*Q1
        Q(4)=Q3+CM3*Q0+CM2*Q1-CM1*Q2
! ---       PROJECTION
        QUAT=1/SQRT(Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2)
        Q(1)=Q(1)*QUAT
        Q(2)=Q(2)*QUAT
        Q(3)=Q(3)*QUAT
        Q(4)=Q(4)*QUAT
! ---       UPDATE M
        POTP = POTENP(Q,PAR)
        CALPHA=HC/CALPHA
        AM1=AM1-H*POTP(1)+FAD1*CM2*CM3*CALPHA
        AM2=AM2-H*POTP(2)+FAD2*CM1*CM3*CALPHA
        AM3=AM3-H*POTP(3)+FAD3*CM1*CM2*CALPHA
      END DO
      AM(1)=AM1+HA*POTP(1)
      AM(2)=AM2+HA*POTP(2)
      AM(3)=AM3+HA*POTP(3)
      RETURN
      
    END SUBROUTINE
    
    PURE SUBROUTINE DMV8(AM,Q,POTENP,H,NSTEP,PAR)
    !--- DISCRETE MOSER-VESELOV ALGORITHM OF ORDER 8
    !--- ORTHOGONAL MATRICES ARE REPRESENTED BY QUATERNIONS
    
      ! Arguments
      procedure(POTENTIAL) :: POTENP
      REAL(rk), intent(inout) :: Q(4),AM(3)
      REAL(rk), intent(in) :: H
      INTEGER, intent(in) :: NSTEP
      type(mmoi), intent(in) :: PAR
    
      ! Subroutine Constants
      integer, PARAMETER :: MSPLIT=1
      
      ! Local Variables
      REAL(rk) :: POTP(3)                                                       ! POTENTIAL ENERGY
      REAL(rk) :: EPS, HA, HC, H2, ERR                                      ! TIME STEP RELATED
      REAL(rk) :: AI1, AI2, AI3, AM1, AM2, AM3                                  ! MMOI VALUES
      REAL(rk) :: XI1, XI2, XI3, XDET, XSI
      REAL(rk) :: XS1, XS2, XS3, XSI12, XSI21
      REAL(rk) :: S3C1, S3C2, D3C1, D3C2, S5C1, S5C2, S5C3,D5C1, D5C2, D5C3
      REAL(rk) :: S7C1, S7C2, S7C3, S7C4, D7C1, D7C2, D7C3, D7C4
      REAL(rk) :: HAM0, HAM2, HAM3
      REAL(rk) :: ANOR0, ANOR2, ANOR3
      REAL(rk) :: ANORHAM, ANOR2HAM, ANORHAM2
      REAL(rk) :: CSS, CDD, AI1MODI, AI2MODI, AI3MODI, AI1MOD, AI2MOD, AI3MOD
      REAL(rk) :: FAD1, FAD2, FAD3, FAC1, FAC2, FAC3
      REAL(rk) :: AM1I, AM2I, AM3I, CM1, CM2, CM3, CM1B, CM2B, CM3B, CALPHA
      REAL(rk) :: Q0, Q1, Q2, Q3, QUAT
      
      INTEGER :: ISTEP, I
      
      EPS=ABS(tiny*H)
      HA=H/2
      HC=4/H
      H2 = H**2
      
      AI1=PAR%AI1
      AI2=PAR%AI2
      AI3=PAR%AI3 
!---      CONSTANT COEFFICIENTS FOR THE MODIFIED MOMENTS OF INERTIA
      XI1=1/AI1+1/AI2+1/AI3
      XI2=1/AI1**2+1/AI2**2+1/AI3**2
      XI3=1/AI1**3+1/AI2**3+1/AI3**3
      XDET=AI1*AI2*AI3
      XS1=AI1+AI2+AI3
      XS2=AI1**2+AI2**2+AI3**2
      XS3=AI1**3+AI2**3+AI3**3
      XSI=(AI1+AI2)/AI3+(AI2+AI3)/AI1+(AI3+AI1)/AI2
      XSI12=(AI1+AI2)/AI3**2+(AI2+AI3)/AI1**2+(AI3+AI1)/AI2**2
      XSI21=(AI1**2+AI2**2)/AI3+(AI2**2+AI3**2)/AI1+(AI3**2+AI1**2)/AI2

      S3C1=XS1/(6*XDET)
      S3C2=-XI1/3
      D3C1=-1/(3*XDET)
      D3C2=S3C1

      S5C1=(XS2/XDET-XI1)/(30*XDET)
      S5C2=(1-XSI)/(30*XDET)
      S5C3=(3*XS1/XDET+2*XI2)/60
      D5C1=-XS1/XDET**2/60
      D5C2=XI1/XDET/10-XS2/XDET**2/60
      D5C3=-(9+XSI)/XDET/60

      S7C1=(4*XDET+17*XS3-15*XDET*XSI)/(2520*XDET**3)
      S7C2=(9*XS1+10*XDET*XI2-6*XSI21)/(420*XDET**2)
      S7C3=((6*XSI12-1.0D2*XI1)*XDET+53*XS2)/(2520*XDET**2)
      S7C4=(15-XDET*XI3-2*XSI)/(630*XDET)
      D7C1=(34*XDET*XI1-19*XS2)/(2520*XDET**3)
      D7C2=(XS3+2*XDET*XSI-85*XDET)/(1260*XDET**3)
      D7C3=(47*XS1+13*XSI21-38*XDET*XI2)/(2520*XDET**2)
      D7C4=(9*XDET*XI1+XDET*XSI12-11*XS2)/(1260*XDET**2)

      POTP = POTENP(Q,PAR)
      AM1=AM(1)-HA*POTP(1)
      AM2=AM(2)-HA*POTP(2)
      AM3=AM(3)-HA*POTP(3)

      DO ISTEP=1,NSTEP
! ---      COMPUTATION OF THE MODIFIED MOMENTS OF INERTIA
        HAM0=0.5_rk*H2*(AM1**2/AI1+AM2**2/AI2+AM3**2/AI3)
        ANOR0=0.5_rk*H2*(AM1**2+AM2**2+AM3**2)
        HAM2=HAM0**2
        ANOR2=ANOR0**2
        HAM3=HAM0**3
        ANOR3=ANOR0**3
        ANORHAM=HAM0*ANOR0
        ANOR2HAM=HAM0*ANOR2
        ANORHAM2=HAM2*ANOR0
        CSS=1+(S3C1*ANOR0+S3C2*HAM0) &
           +(S5C1*ANOR2+S5C2*ANORHAM+S5C3*HAM2) &
           +(S7C1*ANOR3+S7C2*ANOR2HAM+S7C3*ANORHAM2+S7C4*HAM3)
        CDD=(D3C1*ANOR0+D3C2*HAM0) & 
           +(D5C1*ANOR2+D5C2*ANORHAM+D5C3*HAM2) &
           +(D7C1*ANOR3+D7C2*ANOR2HAM+D7C3*ANORHAM2+D7C4*HAM3)
        AI1MODI=CSS/AI1+CDD
        AI2MODI=CSS/AI2+CDD
        AI3MODI=CSS/AI3+CDD
        AI1MOD=1/AI1MODI
        AI2MOD=1/AI2MODI
        AI3MOD=1/AI3MODI
        FAD1=AI2MOD-AI3MOD
        FAD2=AI3MOD-AI1MOD
        FAD3=AI1MOD-AI2MOD
        FAC1=FAD1*AI1MODI
        FAC2=FAD2*AI2MODI
        FAC3=FAD3*AI3MODI
! ---      SOLVE FOR INTERNAL STAGE
        AM1I=AM1*HA*AI1MODI
        AM2I=AM2*HA*AI2MODI
        AM3I=AM3*HA*AI3MODI
        CM1=AM1I+FAC1*AM2I*AM3I
        CM2=AM2I+FAC2*CM1*AM3I
        CM3=AM3I+FAC3*CM1*CM2
        DO I=1,50
          CM1B=CM1
          CM2B=CM2
          CM3B=CM3
          CALPHA=1+CM1**2+CM2**2+CM3**2
          CM1=CALPHA*AM1I+FAC1*CM2*CM3
          CM2=CALPHA*AM2I+FAC2*CM1*CM3
          CM3=CALPHA*AM3I+FAC3*CM1*CM2
          ERR=ABS(CM1B-CM1)+ABS(CM2B-CM2)+ABS(CM3B-CM3)
          IF( ERR < EPS ) EXIT
        END DO
! ---      UPDATE Q
        Q0=Q(1)
        Q1=Q(2)
        Q2=Q(3)
        Q3=Q(4)
        Q(1)=Q0-CM1*Q1-CM2*Q2-CM3*Q3
        Q(2)=Q1+CM1*Q0+CM3*Q2-CM2*Q3
        Q(3)=Q2+CM2*Q0+CM1*Q3-CM3*Q1
        Q(4)=Q3+CM3*Q0+CM2*Q1-CM1*Q2
! ---       PROJECTION
        QUAT=1/SQRT(Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2)
        Q(1)=Q(1)*QUAT
        Q(2)=Q(2)*QUAT
        Q(3)=Q(3)*QUAT
        Q(4)=Q(4)*QUAT
! ---       UPDATE M
        POTP = POTENP(Q,PAR)
        CALPHA=HC/CALPHA
        AM1=AM1-H*POTP(1)+FAD1*CM2*CM3*CALPHA
        AM2=AM2-H*POTP(2)+FAD2*CM1*CM3*CALPHA
        AM3=AM3-H*POTP(3)+FAD3*CM1*CM2*CALPHA
      END DO
      AM(1)=AM1+HA*POTP(1)
      AM(2)=AM2+HA*POTP(2)
      AM(3)=AM3+HA*POTP(3)
      RETURN
            
    END SUBROUTINE    
    
    PURE SUBROUTINE DMV10(AM,Q,POTENP,H,NSTEP,PAR)
    !--- DISCRETE MOSER-VESELOV ALGORITHM OF ORDER 10
    !--- ORTHOGONAL MATRICES ARE REPRESENTED BY QUATERNIONS
        
      ! Arguments
      procedure(POTENTIAL) :: POTENP
      REAL(rk), intent(inout) :: Q(4),AM(3)
      REAL(rk), intent(in) :: H
      INTEGER, intent(in) :: NSTEP
      type(mmoi), intent(in) :: PAR
          
      ! Local Variables
      REAL(rk) :: POTP(3)                                                       ! POTENTIAL ENERGY
      REAL(rk) :: EPS, HA, HC, H2, ERR                                      ! TIME STEP RELATED
      REAL(rk) :: AI1, AI2, AI3, AM1, AM2, AM3                                  ! MMOI VALUES
      REAL(rk) :: XI1, XI2, XI3, XI4, XDET, XDET2, XDET3, XDET4, XSI
      REAL(rk) :: XS1, XS2, XS3, XS4, XSI12, XSI21, XSI13, XSI31, XSI22
      REAL(rk) :: S3C1, S3C2, D3C1, D3C2, S5C1, S5C2, S5C3,D5C1, D5C2, D5C3
      REAL(rk) :: S7C1, S7C2, S7C3, S7C4, D7C1, D7C2, D7C3, D7C4
      REAL(rk) :: S9C1, S9C2 ,S9C3, S9C4, S9C5, D9C1, D9C2, D9C3, D9C4, D9C5      
      REAL(rk) :: HAM0, HAM2, HAM3, HAM4
      REAL(rk) :: ANOR0, ANOR2, ANOR3, ANOR4
      REAL(rk) :: ANORHAM, ANOR2HAM, ANORHAM2, ANOR3HAM, ANOR2HAM2, ANORHAM3
      REAL(rk) :: CSS, CDD, AI1MODI, AI2MODI, AI3MODI, AI1MOD, AI2MOD, AI3MOD
      REAL(rk) :: FAD1, FAD2, FAD3, FAC1, FAC2, FAC3
      REAL(rk) :: AM1I, AM2I, AM3I, CM1, CM2, CM3, CM1B, CM2B, CM3B, CALPHA
      REAL(rk) :: Q0, Q1, Q2, Q3, QUAT
      
      INTEGER :: ISTEP, I
    
      EPS=ABS(tiny*H)
      HA=H/2
      HC=4/H
      H2 = H**2
      AI1=PAR%AI1
      AI2=PAR%AI2
      AI3=PAR%AI3 
    
!--- CONSTANT COEFFICIENTS FOR THE MODIFIED MOMENTS OF INERTIA
      XI1=1/AI1   +1/AI2   +1/AI3
      XI2=1/AI1**2+1/AI2**2+1/AI3**2
      XI3=1/AI1**3+1/AI2**3+1/AI3**3
      XI4=1/AI1**4+1/AI2**4+1/AI3**4
      XDET=AI1*AI2*AI3
      XDET2=XDET**2
      XDET3=XDET**3
      XDET4=XDET2**2
      XS1=AI1+AI2+AI3
      XS2=AI1**2+AI2**2+AI3**2
      XS3=AI1**3+AI2**3+AI3**3
      XS4=AI1**4+AI2**4+AI3**4
      XSI=(AI1+AI2)/AI3+(AI2+AI3)/AI1+(AI3+AI1)/AI2
      XSI12=(AI1+AI2)/AI3**2+(AI2+AI3)/AI1**2+(AI3+AI1)/AI2**2
      XSI21=(AI1**2+AI2**2)/AI3+(AI2**2+AI3**2)/AI1+(AI3**2+AI1**2)/AI2
      XSI13=(AI1+AI2)/AI3**3+(AI2+AI3)/AI1**3+(AI3+AI1)/AI2**3
      XSI31=(AI1**3+AI2**3)/AI3+(AI2**3+AI3**3)/AI1+(AI3**3+AI1**3)/AI2
      XSI22=(AI1**2+AI2**2)/AI3**2+(AI2**2+AI3**2)/AI1**2+(AI3**2+AI1**2)/AI2**2
      
      S3C1=XS1/(6*XDET)
      S3C2=-XI1/3
      D3C1=-1/(3*XDET)
      D3C2=S3C1

      S5C1=(XS2/XDET-XI1)/(30*XDET)
      S5C2=(1-XSI)/(30*XDET)
      S5C3=(3*XS1/XDET+2*XI2)/60
      D5C1=-XS1/XDET2/60
      D5C2=XI1/XDET/10-XS2/XDET2/60
      D5C3=-(9+XSI)/XDET/60

      S7C1=(4*XDET+17*XS3-15*XDET*XSI)/(2520*XDET3)
      S7C2=(9*XS1+10*XDET*XI2-6*XSI21)/(420*XDET2)
      S7C3=((6*XSI12-1.0D2*XI1)*XDET+53*XS2)/(2520*XDET2)
      S7C4=(15-XDET*XI3-2*XSI)/(630*XDET)
      D7C1=(34*XDET*XI1-19*XS2)/(2520*XDET3)
      D7C2=(XS3+2*XDET*XSI-85*XDET)/(1260*XDET3)
      D7C3=(47*XS1+13*XSI21-38*XDET*XI2)/(2520*XDET2)
      D7C4=(9*XDET*XI1+XDET*XSI12-11*XS2)/(1260*XDET2)
      
      S9C1=(62*XS4-94*XDET*XSI21+66*XDET2*XI2+81*XDET*XS1)/(45360*XDET4)
      S9C2=(-77*XSI31+75*XDET*XSI12+214*XS2-240*XDET*XI1)/(22680*XDET3)
      S9C3=(26*XDET*XSI22+55*XS3+204*XDET-50*XDET2*XI3-59*XDET*XSI)/(7560*XDET3)
      S9C4=(137*XDET*XI2-XDET*XSI13+3*XS1-69*XSI21)/(11340*XDET2)
      S9C5=(2*XDET2*XI4+5*XDET*XSI12-171*XDET*XI1+159*XS2)/(45360*XDET2)
      D9C1=(60*XSI*XDET-61*XS3-247*XDET)/(45360*XDET4)
      D9C2=(54*XDET*XS1-XS4+218*XDET*XSI21-426*XDET2*XI2)/(45360*XDET4)
      D9C3=(125*XDET*XI1-5*XSI31-130*XS2+4*XDET*XSI12)/(7560*XDET3)
      D9C4=(67*XS3-735*XDET-15*XDET*XSI22+87*XDET*XSI+34*XDET2*XI3)/(22680*XDET3)
      D9C5=(165*XS1-XDET*XSI13-9*XSI21-145*XDET*XI2)/(45360*XDET2)
      
      POTP = POTENP(Q,PAR)
      
      AM1=AM(1)-HA*POTP(1)
      AM2=AM(2)-HA*POTP(2)
      AM3=AM(3)-HA*POTP(3)
      
      
      DO ISTEP=1,NSTEP
!--- COMPUTATION OF THE MODIFIED MOMENTS OF INERTIA
        HAM0 =0.5_rk*H2*(AM1**2/AI1+AM2**2/AI2+AM3**2/AI3)
        ANOR0=0.5_rk*H2*(AM1**2+AM2**2+AM3**2)
        HAM2=HAM0**2
        ANOR2=ANOR0**2
        HAM3=HAM0**3
        ANOR3=ANOR0**3
        HAM4=HAM2**2
        ANOR4=ANOR2**2
        ANORHAM=HAM0*ANOR0
        ANOR2HAM=HAM0*ANOR2
        ANORHAM2=HAM2*ANOR0
        ANOR3HAM=HAM0*ANOR3
        ANOR2HAM2=HAM2*ANOR2
        ANORHAM3=HAM3*ANOR0
        CSS=1.0_rk+(S3C1*ANOR0+S3C2*HAM0) &
           +(S5C1*ANOR2+S5C2*ANORHAM+S5C3*HAM2) &
           +(S7C1*ANOR3+S7C2*ANOR2HAM+S7C3*ANORHAM2+S7C4*HAM3) &
           +(S9C1*ANOR4+S9C2*ANOR3HAM+S9C3*ANOR2HAM2 &
           +S9C4*ANORHAM3+S9C5*HAM4) 
        CDD=(D3C1*ANOR0+D3C2*HAM0) &
           +(D5C1*ANOR2+D5C2*ANORHAM+D5C3*HAM2) &
           +(D7C1*ANOR3+D7C2*ANOR2HAM+D7C3*ANORHAM2+D7C4*HAM3) &
           +(D9C1*ANOR4+D9C2*ANOR3HAM+D9C3*ANOR2HAM2 &
           +D9C4*ANORHAM3+D9C5*HAM4)
        AI1MODI=CSS/AI1+CDD
        AI2MODI=CSS/AI2+CDD
        AI3MODI=CSS/AI3+CDD
        AI1MOD=1/AI1MODI
        AI2MOD=1/AI2MODI
        AI3MOD=1/AI3MODI
        FAD1=AI2MOD-AI3MOD
        FAD2=AI3MOD-AI1MOD
        FAD3=AI1MOD-AI2MOD
        FAC1=FAD1*AI1MODI
        FAC2=FAD2*AI2MODI
        FAC3=FAD3*AI3MODI        
        
!--- SOLVE FOR INTERNAL STAGE            
        AM1I=AM1*HA*AI1MODI
        AM2I=AM2*HA*AI2MODI
        AM3I=AM3*HA*AI3MODI
        CM1=AM1I+FAC1*AM2I*AM3I
        CM2=AM2I+FAC2*CM1*AM3I
        CM3=AM3I+FAC3*CM1*CM2
        DO I=1, 50
            CM1B=CM1
            CM2B=CM2
            CM3B=CM3
            CALPHA=1.0_rk+CM1**2+CM2**2+CM3**2
            CM1=CALPHA*AM1I+FAC1*CM2*CM3
            CM2=CALPHA*AM2I+FAC2*CM1*CM3
            CM3=CALPHA*AM3I+FAC3*CM1*CM2
            ERR=ABS(CM1B-CM1)+ABS(CM2B-CM2)+ABS(CM3B-CM3)
            IF( ERR < EPS ) EXIT
        END DO 
!--- UPDATE Q
        Q0=Q(1)
        Q1=Q(2)
        Q2=Q(3)
        Q3=Q(4)
        Q(1)=Q0-CM1*Q1-CM2*Q2-CM3*Q3
        Q(2)=Q1+CM1*Q0+CM3*Q2-CM2*Q3
        Q(3)=Q2+CM2*Q0+CM1*Q3-CM3*Q1
        Q(4)=Q3+CM3*Q0+CM2*Q1-CM1*Q2
!--- UPDATE M
        CALPHA=HC/CALPHA
        AM1=AM1+FAD1*CM2*CM3*CALPHA
        AM2=AM2+FAD2*CM1*CM3*CALPHA
        AM3=AM3+FAD3*CM1*CM2*CALPHA                                
        
!--- PROJECTION
        QUAT=1/SQRT(Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2)
        Q(1)=Q(1)*QUAT
        Q(2)=Q(2)*QUAT
        Q(3)=Q(3)*QUAT
        Q(4)=Q(4)*QUAT

        POTP = POTENP(Q,PAR)
        AM1=AM1-H*POTP(1)
        AM2=AM2-H*POTP(2)
        AM3=AM3-H*POTP(3)   
                
      END DO !ISTEP
      AM(1)=AM1+HA*POTP(1)
      AM(2)=AM2+HA*POTP(2)
      AM(3)=AM3+HA*POTP(3)
      
    RETURN      
    END SUBROUTINE

    end module