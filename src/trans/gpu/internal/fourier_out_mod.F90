! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FOURIER_OUT_MOD
CONTAINS
SUBROUTINE FOURIER_OUT(KFIELDS)

!**** *FOURIER_OUT* - Copy fourier data from local array to buffer

!     Purpose.
!     --------
!        Routine for copying fourier data from local array to buffer

!**   Interface.
!     ----------
!     CALL FOURIER_OUT(...)

!     Explicit arguments :  PREEL - local fourier/GP array
!     --------------------  KFIELDS - number of fields
!
!     Externals.  None.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 2000-04-01

!     ------------------------------------------------------------------

USE PARKIND_ECTRANS,ONLY : JPIM     ,JPRBT
USE TPM_DISTR,ONLY : D,MYSETW,MYPROC,NPROC,D_NPTRLS,D_NSTAGTF,D_MSTABF,D_NSTAGT1B,&
	D_NPNTGTB0,D_NPROCM,D_NPROCL
USE TPM_TRANS,ONLY : FOUBUF_IN, ZGTF
USE TPM_GEOMETRY,ONLY : G,G_NMEN,G_NMEN_MAX

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
INTEGER(KIND=JPIM) :: KGL
INTEGER(KIND=JPIM) :: JM,JF,IGLG,IPROC,ISTA
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC, IOFF

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

! like: foubuf_in(1:2*kfields,ista) = zgtf(1:2*kfields,ioff+2*jm)

!$ACC DATA PRESENT(FOUBUF_IN,ZGTF,D_NPTRLS,G_NMEN,D_NPROCM,D_NSTAGT1B,D_MSTABF,D_NPNTGTB0,D_NSTAGTF)&
!$ACC   COPYIN(IBEG,IEND,IINC)
!$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(3) PRIVATE(IGLG,IPROC,ISTA,IOFF)
DO KGL=IBEG,IEND,IINC
   DO JM=0,G_NMEN_MAX      
      DO JF=1,KFIELDS
         IGLG = D_NPTRLS(MYSETW)+KGL-1
         if  (JM <= G_NMEN(IGLG)) then
            IPROC = D_NPROCM(JM)
            ISTA  = (D_NSTAGT1B(D_MSTABF(IPROC))+D_NPNTGTB0(JM,KGL))*2*KFIELDS
            IOFF  = 1+D_NSTAGTF(KGL)
           
            ! imaginary may be not JM+1 but JM+G_NMEN(IGLG)+1
            FOUBUF_IN(ISTA+2*JF-1) = ZGTF(2*JF-1, 2*JM+IOFF)
            FOUBUF_IN(ISTA+2*JF  ) = ZGTF(2*JF  , 2*JM+IOFF)
         end if         
      ENDDO
   ENDDO   
END DO
!$ACC END DATA

END SUBROUTINE FOURIER_OUT
END MODULE FOURIER_OUT_MOD

