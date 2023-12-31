! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE VDTUV_MOD
CONTAINS
SUBROUTINE VDTUV(KFIELD,PEPSNM,PVOR,PDIV,PU,PV)

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,  JPRBT

USE TPM_DIM         ,ONLY : R
USE TPM_FIELDS      ,ONLY : F
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS
use tpm_gen, only: nout


!**** *VDTUV* - Compute U,V in  spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the winds
!        from vorticity and divergence.

!**   Interface.
!     ----------
!        CALL VDTUV(...)

!        Explicit arguments :  KM -zonal wavenumber (input-c)
!        --------------------  KFIELD - number of fields (input-c)
!                              PEPSNM - REPSNM for wavenumber KM (input-c)
!                              PVOR(NLEI1,2*KFIELD) - vorticity (input)
!                              PDIV(NLEI1,2*KFIELD) - divergence (input)
!                              PU(NLEI1,2*KFIELD)   - u wind (output)
!                              PV(NLEI1,2*KFIELD)   - v wind (output)
!        Organisation within NLEI1:
!        NLEI1 = NSMAX+4+mod(NSMAX+4+1,2)
!                        overdimensioning
!        1        : n=NSMAX+2
!        2        : n=NSMAX+1
!        3        : n=NSMAX
!        .        :
!        .        :
!        NSMAX+3  : n=0
!        NSMAX+4  : n=-1

!        Implicit arguments :  Eigenvalues of inverse Laplace operator
!        --------------------  from YOMLAP

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Temperton, 1991, MWR 119 p1303

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-02-01 From VDTUV in IFS CY22R1

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) :: KM, kmloc
INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
REAL(KIND=JPRBT), INTENT(IN)    :: PEPSNM(1:D%NUMP,0:R%NTMAX+2)
REAL(KIND=JPRB), INTENT(IN)    :: PVOR(:,:,:),PDIV(:,:,:)
REAL(KIND=JPRB), INTENT(OUT)   :: PU  (:,:,:),PV  (:,:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: II, IJ, IR, J, JN, ISMAX,JI

!     LOCAL REAL SCALARS
REAL(KIND=JPRBT) :: ZKM
REAL(KIND=JPRBT) :: ZN(-1:R%NTMAX+4)
REAL(KIND=JPRBT) :: ZLAPIN(-1:R%NSMAX+4)
REAL(KIND=JPRBT) :: ZEPSNM(-1:R%NSMAX+4)

!$ACC DATA                                     &
!$ACC      CREATE (ZEPSNM, ZN, ZLAPIN)         &
!$ACC      COPYIN (D,D%MYMS,F,F%RLAPIN,F%RN)   &
!$ACC      PRESENT(PEPSNM, PVOR, PDIV)         &
!$ACC      PRESENT(PU, PV)


!     ------------------------------------------------------------------

!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

ISMAX = R%NSMAX
DO KMLOC=1,D%NUMP
  ZKM = D%MYMS(KMLOC)
  !$ACC PARALLEL LOOP DEFAULT(NONE)
  DO JN=ZKM-1,ISMAX+2
    IJ = ISMAX+3-JN
    ZN(IJ) = F%RN(JN)
    ZLAPIN(IJ) = F%RLAPIN(JN)
    IF( JN >= 0 ) THEN
        ZEPSNM(IJ) = PEPSNM(KMLOC,JN) 
    ELSE
        ZEPSNM(IJ) = 0
    ENDIF
  ENDDO
  !$ACC KERNELS DEFAULT(NONE)
  ZN(0) = F%RN(ISMAX+3)
  !$ACC END KERNELS

!*       1.1      U AND V (KM=0) .

IF(ZKM == 0) THEN
  !$ACC PARALLEL LOOP DEFAULT(NONE)
  DO J=1,KFIELD
    IR = 2*J-1
    DO JI=2,ISMAX+3
      PU(IR,JI,KMLOC) = +&
      &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PVOR(IR,JI+1,KMLOC)-&
      &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PVOR(IR,JI-1,KMLOC)
      PV(IR,JI,KMLOC) = -&
      &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PDIV(IR,JI+1,KMLOC)+&
      &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PDIV(IR,JI-1,KMLOC)
    ENDDO
  ENDDO
ELSE
!*       1.2      U AND V (KM!=0) .

    !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE)
    DO J=1,KFIELD
      DO JI=2,ISMAX+3-ZKM
        !ZKM = D_MYMS(KMLOC)
        IR = 2*J-1
        II = IR+1
        !IF (ZKM>0 .AND. JI<=ISMAX+3-zKM) THEN
          PU(ir,JI,kmloc) = -ZKM*ZLAPIN(JI)*PDIV(ii,JI,kmloc)+&
          &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PVOR(ir,JI+1,kmloc)-&
          &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PVOR(ir,JI-1,kmloc)
          PU(ii,JI,kmloc) = +ZKM*ZLAPIN(JI)*PDIV(ir,JI,kmloc)+&
          &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PVOR(ii,JI+1,kmloc)-&
          &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PVOR(ii,JI-1,kmloc)
          PV(ir,JI,kmloc) = -ZKM*ZLAPIN(JI)*PVOR(ii,JI,kmloc)-&
          &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PDIV(ir,JI+1,kmloc)+&
          &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PDIV(ir,JI-1,kmloc)
          PV(ii,JI,kmloc) = +ZKM*ZLAPIN(JI)*PVOR(ir,JI,kmloc)-&
          &ZN(JI+1)*ZEPSNM(JI)*ZLAPIN(JI+1)*PDIV(ii,JI+1,kmloc)+&
          &ZN(JI-2)*ZEPSNM(JI-1)*ZLAPIN(JI-1)*PDIV(ii,JI-1,kmloc)
        !ENDIF
      ENDDO
    ENDDO
  ENDIF
ENDDO

!$ACC END DATA
!     ------------------------------------------------------------------

END SUBROUTINE VDTUV
END MODULE VDTUV_MOD

