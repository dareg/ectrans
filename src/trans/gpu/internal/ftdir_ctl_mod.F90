! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!**** *FTDIR_CTL - Direct Fourier transform control

!     Purpose. Control routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!     CALL FTDIR_CTL(..)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     PGP     -  gridpoint array
!     KVSETUV - "B" set in spectral/fourier space for
!                u and v variables
!     KVSETSC - "B" set in spectral/fourier space for
!                scalar variables
!     KPTRGP  -  pointer array to fields in gridpoint space

!     Method.
!     -------

!     Externals.  TRGTOL      - transposition routine
!     ----------  FOURIER_OUT - copy fourier data to Fourier buffer
!                 FTDIR       - fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

MODULE FTDIR_CTL_MOD
  USE PARKIND_ECTRANS  ,ONLY : JPIM,JPRB,JPRBT
  USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
CONTAINS
SUBROUTINE FTDIR_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KVSETUV,KVSETSC,KPTRGP,&
 & KVSETSC3A,KVSETSC3B,KVSETSC2,PGP,PGPUV,PGP3A,PGP3B,PGP2)
  USE TPM_GEN, only: nout
  USE TPM_TRANS       ,ONLY : ZGTF,FOUBUF_IN
  USE TRGTOL_MOD      ,ONLY : TRGTOL,TRGTOL_CUDAAWARE
  USE FTDIR_MOD       ,ONLY : FTDIR
  use ieee_arithmetic

  IMPLICIT NONE

  INTEGER(KIND=JPIM),INTENT(IN) :: KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS
  INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
  INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
  INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPTRGP(:)
  INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
  INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
  INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
  REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP(:,:,:)
  REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGPUV(:,:,:,:)
  REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3A(:,:,:,:)
  REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP3B(:,:,:,:)
  REAL(KIND=JPRB),OPTIONAL    , INTENT(IN) :: PGP2(:,:,:)

  INTEGER(KIND=JPIM) :: j3,IVSETUV(KF_UV_G),IVSETSC(KF_SCALARS_G),IVSET(KF_GP),ioff

  IF (PRESENT(KVSETUV)) THEN
	 if (size(KVSETUV) /= KF_UV_G) call abort_trans("Error: kvsetuv")
	 IVSETUV(:) = KVSETUV(:)
  ELSE
	 IVSETUV(:) = -1
  ENDIF

  IF (PRESENT(KVSETSC)) THEN
	 if (size(KVSETSC) /= KF_SCALARS_G) call abort_trans("Error: kvsetsc")
	 IVSETSC(:) = KVSETSC(:)
  ELSE
    IVSETSC(:) = -1
    call initvset(ivsetsc,kf_gp,KF_SCALARS_G,KVSETSC2,KVSETSC3A,KVSETSC3B,PGP3A,PGP3B)
  ENDIF

  ioff = 0

  IF (KF_UV_G > 0) THEN
	 if (2*KF_UV_G > kf_gp) call abor1("Error: ivsetuv")
	 IVSET(ioff+1:ioff+KF_UV_G) = IVSETUV(:)
	 ioff = ioff+KF_UV_G
	 IVSET(ioff+1:ioff+KF_UV_G) = IVSETUV(:)
	 ioff = ioff+KF_UV_G
  ENDIF

  IF (KF_SCALARS_G > 0) THEN
	 if (ioff+kf_scalars_g > kf_gp) call abor1("Error: ivsetuv")
	 IVSET(ioff+1:ioff+KF_SCALARS_G) = IVSETSC(:)
	 ioff = ioff+KF_SCALARS_G
  ENDIF

  write(0,*) "ftdir ivset:",ivset(:)

  ! needed ??? JF_FS=KF_FS-D%IADJUST_D

  !$ACC KERNELS
  ZGTF(:,:) = 0
  !$ACC END KERNELS

#ifdef USE_CUDA_AWARE_MPI_FT
  write(nout,*) "Transpose GP to lat. array (CUDAAWARE)"
  CALL TRGTOL_CUDAAWARE(ZGTF,KF_FS,KF_GP,KF_SCALARS_G,IVSET,KPTRGP,PGP,PGPUV,PGP3A,PGP3B,PGP2)

  ! debug:
  !!$acc update host(zgtf)
#else
  write(nout,*) "Transpose GP to lat. array"
  CALL TRGTOL(ZGTF,KF_FS,KF_GP,KF_SCALARS_G,IVSET,KPTRGP,PGP,PGPUV,PGP3A,PGP3B,PGP2)

  !$ACC UPDATE DEVICE(ZGTF)
#endif

  !do j3=1,size(zgtf,1)
  !  write(nout,*) "zgtf:",minval(zgtf(j3,:)),maxval(zgtf(j3,:))
  !end do

  ! note: no end data (buffer to be passed to ltdir)
  !$ACC ENTER DATA CREATE(FOUBUF_IN)

  write(nout,*) "Direct Fourier transforms - fields:",kf_fs
  if (kf_fs > 0) CALL FTDIR(KF_FS)
END SUBROUTINE FTDIR_CTL

subroutine initvset(ivsetsc,kf_gp,KF_SCALARS_G,KVSETSC2,KVSETSC3A,KVSETSC3B,&
  PGP3A,PGP3B)
  INTEGER(KIND=JPIM),intent(out) :: IVSETSC(:)
  INTEGER(KIND=JPIM),intent(in) :: kf_gp,KF_SCALARS_G
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KVSETSC2(:),KVSETSC3A(:),KVSETSC3B(:)
  REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PGP3A(:,:,:,:)
  REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PGP3B(:,:,:,:)

  INTEGER(KIND=JPIM) :: IFGP2,IFGP3A,IFGP3B,IOFF,J3

  IOFF=0

  IF (PRESENT(KVSETSC2)) THEN
    IFGP2=UBOUND(KVSETSC2,1)
    if (ifgp2 > kf_scalars_g) call abort_trans("Error: kvsetsc2")

    IVSETSC(1:IFGP2)=KVSETSC2(:)
    IOFF=IOFF+IFGP2
  ENDIF

  IF (PRESENT(KVSETSC3A)) THEN
    IFGP3A=UBOUND(KVSETSC3A,1)
    if (ioff+ifgp3a*size(pgp3a,3) > kf_scalars_g) call abort_trans("Error: kvsetsc3a")

    DO J3=1,UBOUND(PGP3A,3)
      IVSETSC(IOFF+1:IOFF+IFGP3A)=KVSETSC3A(:)
      IOFF=IOFF+IFGP3A
    ENDDO
  ENDIF

  IF (PRESENT(KVSETSC3B)) THEN
    IFGP3B=UBOUND(KVSETSC3B,1)
    if (ioff+ifgp3b*size(pgp3b,3) > kf_scalars_g) call abort_trans("Error: kvsetsc3b")

    DO J3=1,UBOUND(PGP3B,3)
      IVSETSC(IOFF+1:IOFF+IFGP3B)=KVSETSC3B(:)
      IOFF=IOFF+IFGP3B
    ENDDO
  ENDIF

  if (ioff /= kf_scalars_g) call abort_trans("Error: ioff /= nf_scalars_g")
end subroutine
END MODULE FTDIR_CTL_MOD
