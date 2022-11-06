# RB_DiscreteMoserVeselov

## Attribution

Original source from https://www.unige.ch/~hairer/software.html

Rigid body integrator, which is described in E. Hairer and G. Vilmart (2006): Preprocessed Discrete Moser-Veselov algorithm for the full dynamics of the rigid body, J. Phys. A: Math. Gen. 39 (2006) 13225-13235.

 - `DMV10.f` - Source file for DMV solvers of order 2,4,6,8, and 10.
 - `DR_DMV10.f` - Driver program for DMV solver with example code simulating for 10 seconds in 1000 steps.
 
 ## Derived Works
 
  - `DMV_F77/DMV10.FI` - Original Fortrran 77 solver subroutines as an Fortran include file.
  - `DMV_F77/DMV_DRIVER.F` - Original Fortrran 77 example code as an Fortran source file.
  - `DMV_F90/mod_dmv.f90` - Ported Fortran 95 solver subroutines in a Fortran module
  - `DMV_F90/DMV_DRIVER.f90` - Ported Fortran 95 example code as a Fortran program.
  
 # Source Code
 
 Here is the description from the original Fortran 77 files.

```fortran
C-----------------------------------------------------------------------
C            PREPROCESSED DISCRETE MOSER-VESELOV ALGORITHM 
C-----------------------------------------------------------------------
C  PREPROCESSED DISCRETE MOSER-VESELOV ALGORITHM OF ORDER 10 FOR THE 
C  NUMERICAL SOLUTION OF THE EQUATIONS OF MOTION OF THE FREE RIGID BODY. 
C  ORTHOGONAL MATRICES ARE REPRESENTED BY QUATERNIONS. 
C  THE CODE IS READY TO INCLUDE AN EXERNAL POTENTIAL (SYMMETRIC STRANG  
C  SPLITTING OF ORDER 2). 
C-----------------------------------------------------------------------
C     AUTHORS: ERNST HAIRER (1) AND GILLES VILMART (1)(2)
C                (1) UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES
C                    2-4 RUE DU LIEVRE, CASE POSTALE 64
C                    CH-1211 GENEVE 4, SWITZERLAND 
C                (2) IRISA/INRIA RENNES, PROJET IPSO
C                    CAMPUS DE BEAULIEU, F-35042 RENNES CEDEX, FRANCE
C     E-MAILS: Ernst.Hairer@math.unige.ch
C              Gilles.Vilmart@math.unige.ch
C
C  THIS CODE IS DESCRIBED IN:
C    E. HAIRER AND G. VILMART, PREPROCESSED DISCRETE MOSER-VESELOV
C    ALGORITHM FOR THE FULL DYNAMICS OF A RIGID BODY
C    J. Phys. A: Math. Gen. 39 (2006) 13225-13235. 
C    http://stacks.iop.org/0305-4470/39/13225
C    https://www.unige.ch/~hairer/preprints/predmv.pdf
C
C     VERSION: AUGUST 29, 2006
C     (latest correction of a small bug: December 28, 2007)
C-----------------------------------------------------------------------
C     INPUT PARAMETERS
C     ----------------
C     AM(I)      INITIAL ANGULAR MOMENTUM (I=1,2,3)
C     Q(I)       INITIAL QUATERNION FOR ORTHOGONAL MATRIX (I=1,2,3,4) 
C     H          STEP SIZE
C     NSTEP      NUMBER OF STEPS
C     POTENP     NAME (EXTERNAL) OF SUBROUTINE FOR AN EXTERNAL POTENTIAL
C                   SUBROUTINE POTENP(Q,POTP,RPAR,IPAR)
C                   DIMENSION Q(4),POTP(3)
C                   POTP(1)=...   ETC.
C     RPAR,IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH
C                CAN BE USED FOR COMMUNICATION BETWEEN SUBROUTINES
C     RPAR(11), RPAR(12), RPAR(13) ARE THE THREE MOMENTS OF INERTIA OF
C                THE RIGID BODY
C
C     OUTPUT PARAMETERS
C     -----------------
C     AM(I)       SOLUTION (ANGULAR MOMENTUM) AT ENDPOINT
C     Q(I)        SOLUTION (QUATERNION) AT ENDPOINT
C----------------------------------------------------------------------- 
```
