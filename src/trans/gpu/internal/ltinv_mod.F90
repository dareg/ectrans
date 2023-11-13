! (C) Copyright 2000- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE LTINV_MOD
  USE PARKIND1  ,ONLY : JPIM     ,JPRB
  USE PARKIND_ECTRANS  ,ONLY : JPRBT
  USE TPM_DIM         ,ONLY : R,R_NSMAX, R_NDGNH, R_NDGL
  USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS,D_NASM0,D_NSTAGT0B,D_NSTAGT1B,D_NPROCL,D_NPNTGTB1
  USE TPM_TRANS       ,ONLY : LDIVGP, LVORGP, NF_SC2, NF_SC3A, NF_SC3B, foubuf_in
  USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
  USE TPM_FIELDS      ,ONLY : F,ZIA,ZSOA1,ZAOA1,ZEPSNM
  USE YOMHOOK,ONLY : LHOOK, DR_HOOK, JPHOOK

  IMPLICIT NONE

CONTAINS
SUBROUTINE LTINV(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,KLEI2,KDIM1,&
 & PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)
  !USE TPM_FLT
  !USE TPM_GEOMETRY
  USE VDTUV_MOD       ,ONLY : VDTUV
  USE SPNSDE_MOD      ,ONLY : SPNSDE
  USE LEINV_MOD       ,ONLY : LEINV
  USE FSPGL_INT_MOD   ,ONLY : FSPGL_INT
  USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

  !**** *LTINV* - Inverse Legendre transform
  !
  !     Purpose.
  !     --------
  !        Tranform from Laplace space to Fourier space, compute U and V
  !        and north/south derivatives of state variables.

  !**   Interface.
  !     ----------
  !        *CALL* *LTINV(...)

  !        Explicit arguments :
  !        --------------------
  !          KM        - zonal wavenumber
  !          KMLOC     - local zonal wavenumber
  !          PSPVOR    - spectral vorticity
  !          PSPDIV    - spectral divergence
  !          PSPSCALAR - spectral scalar variables

  !        Implicit arguments :  The Laplace arrays of the model.
  !        --------------------  The values of the Legendre polynomials
  !                              The grid point arrays of the model
  !     Method.
  !     -------

  !     Externals.
  !     ----------

  !         PREPSNM - prepare REPSNM for wavenumber KM
  !         PRFI1B  - prepares the spectral fields
  !         VDTUV   - compute u and v from vorticity and divergence
  !         SPNSDE  - compute north-south derivatives
  !         LEINV   - Inverse Legendre transform
  !         ASRE1   - recombination of symmetric/antisymmetric part

  !     Reference.
  !     ----------
  !        ECMWF Research Department documentation of the IFS
  !        Temperton, 1991, MWR 119 p1303

  !     Author.
  !     -------
  !        Mats Hamrud  *ECMWF*

  !     Modifications.
  !     --------------
  !        Original : 00-02-01 From LTINV in IFS CY22R1
  !     ------------------------------------------------------------------

  INTEGER(KIND=JPIM), INTENT(IN) :: KF_OUT_LT
  INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
  INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
  INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCDERS
  INTEGER(KIND=JPIM), INTENT(IN) :: KLEI2
  INTEGER(KIND=JPIM), INTENT(IN) :: KDIM1
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPVOR(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPDIV(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSCALAR(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC2(:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC3A(:,:,:)
  REAL(KIND=JPRB)   ,OPTIONAL,INTENT(IN)  :: PSPSC3B(:,:,:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRUV(:)
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)  :: KFLDPTRSC(:)

  EXTERNAL  FSPGL_PROC
  OPTIONAL  FSPGL_PROC

  INTEGER(KIND=JPIM) :: IFC, ISTA, IIFC, IDGLU
  INTEGER(KIND=JPIM) :: IVORL,IVORU,IDIVL,IDIVU,IUL,IUU,IVL,IVU,ISL,ISU,IDL,IDU
  INTEGER(KIND=JPIM) :: IFIRST, ILAST, IDIM1,IDIM2,IDIM3,J3
  INTEGER(KIND=JPIM) :: KM,KMLOC
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  ! added
  IF (PRESENT(KFLDPTRUV)) call abort_trans("Error: kfldptruv not supported in GPU version")
  IF (PRESENT(KFLDPTRSC)) call abort_trans("Error: kfldptrsc not supported in GPU version")

  IF (LHOOK) CALL DR_HOOK('LTINV_MOD',0,ZHOOK_HANDLE)

  !*       3.    SPECTRAL COMPUTATIONS FOR U,V AND DERIVATIVES.
  !              ----------------------------------------------

  !ILEI2 = 2*(4*KF_UV+KF_SCALARS+KF_SCDERS)
  !IDIM1 = 2*KF_OUT_LT

  IFIRST = 1
  ILAST  = 0

  !*       1.    PREPARE ZEPSNM.
  !              ---------------

  IF (KF_UV > 0) THEN
    IVORL = 1
    IVORU = 2*KF_UV
    IDIVL = 2*KF_UV+1
    IDIVU = 4*KF_UV
    IUL   = 4*KF_UV+1
    IUU   = 6*KF_UV
    IVL   = 6*KF_UV+1
    IVU   = 8*KF_UV

    IDIM2=UBOUND(PSPVOR,2)

    !$ACC DATA COPYIN(PSPVOR,PSPDIV)

    !CALL PRFI1B(ZIA,PSPVOR,KF_UV,IDIM2,KFLDPTRUV,ivorl-1)
    !CALL PRFI1B(ZIA,PSPDIV,KF_UV,IDIM2,KFLDPTRUV,idivl-1)
    CALL PRFI1A(ZIA,PSPVOR,KF_UV,IDIM2,ivorl-1)
    CALL PRFI1A(ZIA,PSPDIV,KF_UV,IDIM2,idivl-1)

    !$ACC END DATA

    CALL VDTUV(KF_UV,ZEPSNM,ZIA(IVORL:IVORU,:,:),ZIA(IDIVL:IDIVU,:,:),&
     & ZIA(IUL:IUU,:,:),ZIA(IVL:IVU,:,:))
    ILAST = ILAST+8*KF_UV
  ENDIF

  IF(KF_SCALARS > 0)THEN
    IF(PRESENT(PSPSCALAR)) THEN
      IFIRST = ILAST+1
      ILAST  = IFIRST-1+2*KF_SCALARS

      IDIM2=UBOUND(PSPSCALAR,2)
      !$ACC DATA COPYIN(PSPSCALAR)
      !CALL PRFI1B(ZIA,PSPSCALAR,KF_SCALARS,IDIM2,KFLDPTRSC,ifirst-1)
      CALL PRFI1A(ZIA,PSPSCALAR,KF_SCALARS,IDIM2,ifirst-1)
      !$ACC END DATA
    ELSE
      IF(PRESENT(PSPSC2) .AND. NF_SC2 > 0) THEN
        IFIRST = ILAST+1
        ILAST  = IFIRST-1+2*NF_SC2
        IDIM2=UBOUND(PSPSC2,2)
        !$ACC DATA COPYIN(PSPSC2)
        CALL PRFI1A(ZIA,PSPSC2,NF_SC2,IDIM2,ifirst-1)
        !$ACC END DATA
      ENDIF
      IF(PRESENT(PSPSC3A) .AND. NF_SC3A > 0) THEN
        IDIM1=NF_SC3A
        IDIM3=UBOUND(PSPSC3A,3)
        IDIM2=UBOUND(PSPSC3A,2)
        !$ACC DATA COPYIN(PSPSC3A)
        DO J3=1,IDIM3
          IFIRST = ILAST+1
          ILAST  = IFIRST-1+2*IDIM1

          CALL PRFI1A(ZIA,PSPSC3A(:,:,J3),IDIM1,IDIM2,ifirst-1)
        ENDDO
        !$ACC END DATA
      ENDIF
      IF(PRESENT(PSPSC3B) .AND. NF_SC3B > 0) THEN
        IDIM1=NF_SC3B
        IDIM3=UBOUND(PSPSC3B,3)
        IDIM2=UBOUND(PSPSC3B,2)
        !$ACC DATA COPYIN(PSPSC3B)
        DO J3=1,IDIM3
          IFIRST = ILAST+1
          ILAST  = IFIRST-1+2*IDIM1

          CALL PRFI1A(ZIA,PSPSC3B(:,:,J3),IDIM1,IDIM2,ifirst-1)
        ENDDO
        !$ACC END DATA
      ENDIF
    ENDIF

    IF(ILAST /= 8*KF_UV+2*KF_SCALARS) THEN
      WRITE(0,*) 'LTINV:KF_UV,KF_SCALARS,ILAST ',KF_UV,KF_SCALARS,ILAST
      CALL ABORT_TRANS('LTINV_MOD:ILAST /= 8*KF_UV+2*KF_SCALARS')
    ENDIF
  ENDIF

  IF (KF_SCDERS > 0) THEN
     ! stop 'Error: code path not (yet) supported in GPU version'
     ISL = 2*(4*KF_UV)+1
     ISU = ISL+2*KF_SCALARS-1
     IDL = 2*(4*KF_UV+KF_SCALARS)+1
     IDU = IDL+2*KF_SCDERS-1
     write(0,*) "ltinv ders (spnsde) - inds:",isl,isu,idl,idu
     CALL SPNSDE(KF_SCALARS,ZEPSNM,ZIA(ISL:ISU,:,:),ZIA(IDL:IDU,:,:))
  ENDIF

  !*       4.    INVERSE LEGENDRE TRANSFORM.
  !              ---------------------------

  ISTA = 1
  IFC  = 2*KF_OUT_LT
  IF(KF_UV > 0 .AND. .NOT. LVORGP) ISTA = ISTA+2*KF_UV
  IF(KF_UV > 0 .AND. .NOT. LDIVGP) ISTA = ISTA+2*KF_UV

  IF( KF_OUT_LT > 0 ) THEN
    CALL LEINV(IFC,KF_OUT_LT,ZIA(ISTA:ISTA+IFC-1,:,:),ZAOA1,ZSOA1)

    ! RECOMBINATION SYMMETRIC/ANTISYMMETRIC PART.
    CALL ASRE1B(KF_OUT_LT,ZAOA1,ZSOA1)

    ! OPTIONAL COMPUTATIONS IN FOURIER SPACE
    IF(PRESENT(FSPGL_PROC)) THEN
      stop 'Error: SPGL_PROC is not (yet) optimized in GPU version'
      CALL FSPGL_INT(KF_UV,KF_SCALARS,KF_SCDERS,KF_OUT_LT,FSPGL_PROC,KFLDPTRUV,KFLDPTRSC)
    ENDIF
  ENDIF

  IF (LHOOK) CALL DR_HOOK('LTINV_MOD',1,ZHOOK_HANDLE)
END SUBROUTINE LTINV

SUBROUTINE PRFI1A(PIA,PSPEC,KFIELDS,KDIM,ioff)
  INTEGER(KIND=JPIM),INTENT(IN)   :: KFIELDS
  REAL(KIND=JPRB)   ,INTENT(IN)   :: PSPEC(:,:)
  REAL(KIND=JPRB)   ,INTENT(OUT)  :: PIA(:,:,:)
  INTEGER(KIND=JPIM),INTENT(IN) :: KDIM,ioff

  INTEGER(KIND=JPIM) :: KM,KMLOC,INM, JN, JFLD

  !$ACC DATA copyin(ioff) PRESENT(D_NUMP,R_NSMAX,D_MYMS,D_NASM0,PIA,PSPEC)

  !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,INM)
  DO KMLOC=1,D_NUMP
     DO JN=0,R_NSMAX
        DO JFLD=1,KFIELDS
           KM = D_MYMS(KMLOC)
           if (JN <= R_NSMAX-KM) then
              INM = D_NASM0(KM)+2*(R_NSMAX-KM-JN)
              IF( INM < KDIM ) THEN
              PIA(ioff+2*JFLD-1,JN+3,KMLOC) = PSPEC(JFLD,INM)
              PIA(ioff+2*JFLD,JN+3,KMLOC) = PSPEC(JFLD,INM+1)
              ENDIF
           end if
        ENDDO
     ENDDO
  END DO

  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM)
  DO KMLOC=1,D_NUMP
     DO JFLD=ioff+1,ioff+2*KFIELDS
        KM = D_MYMS(KMLOC)
        PIA(JFLD,1,KMLOC) = 0.0_JPRB
        PIA(JFLD,2,KMLOC) = 0.0_JPRB
        PIA(JFLD,R_NSMAX+4-KM,KMLOC) = 0.0_JPRB
     ENDDO
  END DO

  !$ACC END DATA
END SUBROUTINE

SUBROUTINE PRFI1B(PIA,PSPEC,KFIELDS,KDIM,KFLDPTR,ioff)
   INTEGER(KIND=JPIM),INTENT(IN) :: KFIELDS
   REAL(KIND=JPRB)   ,INTENT(IN)   :: PSPEC(:,:)
   REAL(KIND=JPRB)   ,INTENT(OUT)  :: PIA(:,:,:)
   INTEGER(KIND=JPIM),INTENT(IN) :: KDIM,ioff
   INTEGER(KIND=JPIM),INTENT(IN) :: KFLDPTR(:)

   INTEGER(KIND=JPIM) :: KM,KMLOC,INM, JN, JFLD

   !$ACC DATA copyin(ioff) PRESENT(D_NUMP,R_NSMAX,D_MYMS,D_NASM0,PIA,PSPEC,KFLDPTR)

   !$ACC PARALLEL LOOP COLLAPSE(3) DEFAULT(NONE) PRIVATE(KM,INM)
   DO KMLOC=1,D_NUMP
      DO JN=0,R_NSMAX
         DO JFLD=1,KFIELDS
            KM = D_MYMS(KMLOC)
            IF (JN <= R_NSMAX-KM) THEN
               INM = D_NASM0(KM)+2*(R_NSMAX-KM-JN)
               PIA(ioff+2*JFLD-1,JN+3,KMLOC) = PSPEC(KFLDPTR(JFLD),INM)
               PIA(ioff+2*JFLD,JN+3,KMLOC) = PSPEC(KFLDPTR(JFLD),INM+1)
            END IF
         ENDDO
      ENDDO
   END DO

   !$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(2) PRIVATE(KM)
   DO KMLOC=1,D_NUMP
      DO JFLD=ioff+1,ioff+2*KFIELDS
         KM = D_MYMS(KMLOC)
         PIA(JFLD,1,KMLOC) = 0.0_JPRB
         PIA(JFLD,2,KMLOC) = 0.0_JPRB
         PIA(JFLD,R_NSMAX+4-KM,KMLOC) = 0.0_JPRB
      ENDDO
   END DO

   !$ACC END DATA
END SUBROUTINE

SUBROUTINE ASRE1B(KFIELD,PAOA,PSOA)
   INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
   REAL(KIND=JPRBT),   INTENT(IN)  :: PSOA(:,:,:)
   REAL(KIND=JPRBT),   INTENT(IN)  :: PAOA(:,:,:)

   INTEGER(KIND=JPIM) :: KM,KMLOC,ISL, IGLS, JFLD, JGL , ISTAN, ISTAS

   !$ACC DATA PRESENT(PAOA,PSOA,D_MYMS,D_NPROCL,D_NSTAGT0B,D_NPNTGTB1,G_NDGLU,FOUBUF_IN)

   !$ACC PARALLEL LOOP DEFAULT(NONE) PRIVATE(KM,ISL,JGL,JFLD,ISTAN,IGLS,ISTAS)
   DO KMLOC=1,D_NUMP
      KM = D_MYMS(KMLOC)
      ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)
      DO JGL=ISL, R_NDGNH
         ISTAN = (D_NSTAGT0B(D_NPROCL(JGL)) + D_NPNTGTB1(KMLOC,JGL))*2*KFIELD
         IGLS = R_NDGL+1-JGL
         ISTAS = (D_NSTAGT0B(D_NPROCL(IGLS)) + D_NPNTGTB1(KMLOC,IGLS))*2*KFIELD
         DO JFLD=1,2*KFIELD
            FOUBUF_IN(ISTAN+JFLD) = PSOA(JFLD,JGL,KMLOC)+PAOA(JFLD,JGL,KMLOC)
            FOUBUF_IN(ISTAS+JFLD) = PSOA(JFLD,JGL,KMLOC)-PAOA(JFLD,JGL,KMLOC)
         ENDDO
      ENDDO
   ENDDO

   !$ACC END DATA
END SUBROUTINE
END MODULE LTINV_MOD

