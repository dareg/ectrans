! (C) Copyright 1987- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTDIR_MOD
CONTAINS
SUBROUTINE LTDIR(KF_FS,KF_UV,KF_SCALARS,KLED2,PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2,KFLDPTRUV,KFLDPTRSC)
  USE PARKIND1  ,ONLY : JPIM     ,JPRB
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
  use tpm_gen,only: nout,nprintlev
  USE TPM_DIM     ,ONLY : R
  USE TPM_DISTR   ,ONLY : D
  USE TPM_GEOMETRY
  USE PREPSNM_MOD ,ONLY : PREPSNM
  USE PRFI2B_MOD  ,ONLY : PRFI2B
  USE LDFOU2_MOD  ,ONLY : LDFOU2
  USE LEDIR_MOD   ,ONLY : LEDIR
  USE UVTVD_MOD
  USE UPDSP_MOD   ,ONLY : UPDSP
  USE TPM_FIELDS      ,ONLY : ZAIA,ZOA1,ZOA2,ZEPSNM

  !**** *LTDIR* - Control of Direct Legendre transform step

  !     Purpose.
  !     --------
  !        Tranform from Fourier space to spectral space, compute
  !        vorticity and divergence.

  !**   Interface.
  !     ----------
  !        *CALL* *LTDIR(...)*

  !        Explicit arguments :
  !        --------------------  KM     - zonal wavenumber
  !                              KMLOC  - local zonal wavenumber

  !        Implicit arguments :  None
  !        --------------------

  !     Method.
  !     -------

  !     Externals.
  !     ----------
  !         PREPSNM - prepare REPSNM for wavenumber KM
  !         PRFI2   - prepares the Fourier work arrays for model variables.
  !         LDFOU2  - computations in Fourier space
  !         LEDIR   - direct Legendre transform
  !         UVTVD   -
  !         UPDSP   - updating of spectral arrays (fields)

  !     Reference.
  !     ----------
  !        ECMWF Research Department documentation of the IFS

  !     Author.
  !     -------
  !        Mats Hamrud and Philippe Courtier  *ECMWF*

  !     Modifications.
  !     --------------
  !        Original : 87-11-24
  !        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
  !                            for uv formulation
  !        Modified 93-03-19 D. Giard - CDCONF='T' for tendencies
  !        Modified 93-11-18 M. Hamrud - use only one Fourier buffer
  !        Modified 94-04-06 R. El khatib Full-POS implementation
  !        M.Hamrud  : 94-11-01 New conf 'G' - vor,div->vor,div
  !                             instead of u,v->vor,div
  !        MPP Group : 95-10-01 Support for Distributed Memory version
  !        K. YESSAD (AUGUST 1996):
  !               - Legendre transforms for transmission coefficients.
  !        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
  !        R. El Khatib 12-Jul-2012 LDSPC2 replaced by UVTVD
  !     ------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER(KIND=JPIM),INTENT(IN)   :: KF_FS,KF_UV,KF_SCALARS,KLED2
  REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
  REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
  REAL(KIND=JPRB)  ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC2(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3A(:,:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT) :: PSPSC3B(:,:,:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)

  INTEGER(KIND=JPIM) :: IUS, IUE, IVS, IVE, IVORS, IVORE, IDIVS, IDIVE
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('LTDIR_MOD',0,ZHOOK_HANDLE)

  !*       1.    PREPARE LEGENDRE POLONOMIALS AND EPSNM
  !              --------------------------------------

  write(nout,*) "anti-sym. ledir - nfs/nuv/nle:",kf_fs,kf_uv,kled2
  CALL PRFI2B(KF_FS,ZAIA,-1)
  CALL LDFOU2(KF_UV,ZAIA)
  CALL LEDIR(KF_FS,KLED2,ZAIA,ZOA1,-1)

  write(nout,*) "sym. ledir (idem)"
  CALL PRFI2B(KF_FS,ZAIA,1)
  CALL LDFOU2(KF_UV,ZAIA)
  CALL LEDIR(KF_FS,KLED2,ZAIA,ZOA1,1)

  !*       5.    COMPUTE VORTICITY AND DIVERGENCE.
  !              ---------------------------------

  IF (KF_UV > 0) THEN
     write(nout,*) "call uvtvd" ,kf_uv
     IUS = 1
     IUE = 2*KF_UV
     IVS = 2*KF_UV+1
     IVE = 4*KF_UV
     IVORS = 1
     IVORE = 2*KF_UV
     IDIVS = 2*KF_UV+1
     IDIVE = 4*KF_UV
     CALL UVTVD(KF_UV)
  ENDIF

  !*       6.    UPDATE SPECTRAL ARRAYS.
  !              -----------------------

  ! this is on the host, so need to cp from device, Nils
  if (nprintlev > 0) write(nout,*) "update spec. (updsp)",KF_UV,KF_SCALARS
  CALL UPDSP(KF_UV,KF_SCALARS,ZOA1,ZOA2,PSPVOR,PSPDIV,PSPSCALAR,&
   & PSPSC3A,PSPSC3B,PSPSC2,KFLDPTRUV,KFLDPTRSC)

  IF (LHOOK) CALL DR_HOOK('LTDIR_MOD',1,ZHOOK_HANDLE)
  END SUBROUTINE LTDIR
  END MODULE LTDIR_MOD
