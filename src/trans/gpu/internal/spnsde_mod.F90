! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SPNSDE_MOD
CONTAINS
SUBROUTINE SPNSDE(KF_SCALARS,PEPSNM,PF,PNSD)

USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRB,  JPRBT

USE TPM_GEN, only: nout
USE TPM_DIM         ,ONLY : R
USE TPM_FIELDS      ,ONLY : F
USE TPM_DISTR       ,ONLY : D
!USE TPM_TRANS


!**** *SPNSDE* - Compute North-South derivative in spectral space

!     Purpose.
!     --------
!        In Laplace space compute the the North-south derivative

!**   Interface.
!     ----------
!        CALL SPNSDE(...)

!        Explicit arguments :
!        --------------------
!        KM -zonal wavenumber (input-c)
!        PEPSNM - REPSNM for wavenumber KM (input-c)
!        PF  (NLEI1,2*KF_SCALARS) - input field (input)
!        PNSD(NLEI1,2*KF_SCALARS) - N-S derivative (output)

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

!        Implicit arguments :  YOMLAP
!        --------------------

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
!        Original : 00-02-01 From SPNSDE in IFS CY22R1

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM)  :: KM, KMLOC
INTEGER(KIND=JPIM), INTENT(IN)  :: KF_SCALARS
!REAL(KIND=JPRBT),    INTENT(IN)  :: PEPSNM(0:R%NTMAX+2)
REAL(KIND=JPRBT),    INTENT(IN)  :: PEPSNM(1:D%NUMP,0:R%NTMAX+2)
REAL(KIND=JPRB),    INTENT(IN)  :: PF(:,:,:)
REAL(KIND=JPRB),    INTENT(OUT) :: PNSD(:,:,:)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: IJ, ISKIP, J, JN,JI,ISMAX, IR, II
REAL(KIND=JPRBT) :: ZZEPSNM(-1:R%NSMAX+4)
REAL(KIND=JPRBT) :: ZN(-1:R%NTMAX+4)

!$ACC DATA CREATE (ZN,ZZEPSNM) PRESENT (F,F%RN,PEPSNM,PF,PNSD)

!     ------------------------------------------------------------------

!*       1.    COMPUTE NORTH SOUTH DERIVATIVE.
!              -------------------------------


!*       1.1      COMPUTE

ISMAX = R%NSMAX
!loop over wavenumber
DO KMLOC=1,D%NUMP
  KM = D%MYMS(KMLOC)
  !$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(IJ)
  DO JN=KM-1,ISMAX+2
   IJ = ISMAX+3-JN
   ZN(IJ) = F%RN(JN)
   IF( JN >= 0 ) THEN
       ZZEPSNM(IJ) = PEPSNM(KMLOC,JN)
   ELSE
       ZZEPSNM(IJ) = 0
   ENDIF
  ENDDO

  !$ACC KERNELS DEFAULT(NONE)
  ZN(0) = F%RN(ISMAX+3)
  !$ACC END KERNELS

  IF(KM == 0) THEN
      !$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(IR,ji)
      DO J=1,KF_SCALARS
        IR = 2*J-1
        DO JI=2,ISMAX+3
          PNSD(IR,JI,KMLOC) = -ZN(JI+1)*ZZEPSNM(JI)*PF(IR,JI+1,KMLOC)+&
            &ZN(JI-2)*ZZEPSNM(JI-1)*PF(IR,JI-1,KMLOC)
        ENDDO
      ENDDO
  ELSE  
    !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE) PRIVATE(ji,IR,II)
    DO J=1,KF_SCALARS
      DO JI=2,ISMAX+3-KM
        IR = 2*J-1
        II = IR+1
        PNSD(IR,JI,KMLOC) = -ZN(JI+1)*ZZEPSNM(JI)*PF(IR,JI+1,KMLOC)+&
          &ZN(JI-2)*ZZEPSNM(JI-1)*PF(IR,JI-1,KMLOC)
        PNSD(II,JI,KMLOC) = -ZN(JI+1)*ZZEPSNM(JI)*PF(II,JI+1,KMLOC)+&
          &ZN(JI-2)*ZZEPSNM(JI-1)*PF(II,JI-1,KMLOC)
      ENDDO
    ENDDO
  ENDIF
END DO

!$ACC END DATA
END SUBROUTINE SPNSDE
END MODULE SPNSDE_MOD
