! (C) Copyright 2001- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE DIR_TRANS_CTL_MOD
USE PARKIND1  ,ONLY : JPIM     ,JPRB

CONTAINS
SUBROUTINE DIR_TRANS_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KF_UV,KF_SCALARS,&
 & PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,&
 & PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2)

!**** *DIR_TRANS_CTL* - Control routine for direct spectral transform.

!     Purpose.
!     --------
!        Control routine for the direct spectral transform

!**   Interface.
!     ----------
!     CALL DIR_TRANS_CTL(...)

!     Explicit arguments :
!     --------------------
!     KF_UV_G      - global number of spectral u-v fields
!     KF_SCALARS_G - global number of scalar spectral fields
!     KF_GP        - total number of output gridpoint fields
!     KF_FS        - total number of fields in fourier space
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity
!     PSPDIV(:,:)  - spectral divergence
!     PSPSCALAR(:,:) - spectral scalarvalued fields
!     KVSETUV(:)  - indicating which 'b-set' in spectral space owns a
!                   vor/div field. Equivalant to NBSETLEV in the IFS.
!                   The length of KVSETUV should be the GLOBAL number
!                   of u/v fields which is the dimension of u and v releated
!                   fields in grid-point space.
!     KVESETSC(:) - indicating which 'b-set' in spectral space owns a
!                   scalar field. As for KVSETUV this argument is required
!                   if the total number of processors is greater than
!                   the number of processors used for distribution in
!                   spectral wave space.
!     PGP(:,:,:)  - gridpoint fields

!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):
!
!       u             : KF_UV_G fields
!       v             : KF_UV_G fields
!       scalar fields : KF_SCALARS_G fields

!     Method.
!     -------

!     Externals.  SHUFFLE     - reshuffle fields for load balancing
!     ----------  FIELD_SPLIT - split fields in NPROMATR packets
!                 LTDIR_CTL   - control of Legendre transform
!                 FTDIR_CTL   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-01-03

!     ------------------------------------------------------------------

USE TPM_GEN         ,ONLY : nout,NPROMATR
USE TPM_TRANS       ,ONLY : FOUBUF_IN, NF_SC2, NF_SC3A, NF_SC3B
USE LTDIR_CTL_MOD   ,ONLY : LTDIR_CTL
USE FTDIR_CTL_MOD   ,ONLY : FTDIR_CTL
USE SHUFFLE_MOD     ,ONLY : SHUFFLE
USE FIELD_SPLIT_MOD ,ONLY : FIELD_SPLIT

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS_G
INTEGER(KIND=JPIM), INTENT(IN) :: KF_GP
INTEGER(KIND=JPIM), INTENT(IN) :: KF_FS
INTEGER(KIND=JPIM), INTENT(IN) :: KF_UV
INTEGER(KIND=JPIM), INTENT(IN) :: KF_SCALARS
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSETSC2(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGPUV(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP3A(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP3B(:,:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PGP2(:,:,:)

INTEGER(KIND=JPIM) :: IPTRGP(KF_GP),IPTRSPUV(NPROMATR),IPTRSPSC(NPROMATR)
INTEGER(KIND=JPIM) :: ISHFUV_G(KF_GP),ISHFSC_G(KF_GP)
INTEGER(KIND=JPIM) :: IVSETUV(KF_GP),IVSETSC(KF_GP)
INTEGER(KIND=JPIM) :: IBLKS,JBLK,ISTUV_G,IENUV_G
INTEGER(KIND=JPIM) :: IF_UV_G,IF_UV,ISTUV,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER(KIND=JPIM) :: jf,JFLD,ISTSC_G,IENSC_G,ISTSC,IENSC,IENUV,IF_GPB

! Perform transform

IF_GPB = 2*KF_UV_G+KF_SCALARS_G

IF(NPROMATR > 0 .AND. IF_GPB > NPROMATR) THEN
  ! Fields to be split into packets

  CALL SHUFFLE(KF_UV_G,KF_SCALARS_G,ISHFUV_G,IVSETUV,ISHFSC_G,IVSETSC,KVSETUV,KVSETSC)

  IBLKS=(IF_GPB-1)/NPROMATR+1

  DO JBLK=1,IBLKS
    CALL FIELD_SPLIT(JBLK,KF_GP,KF_UV_G,IVSETUV,IVSETSC,&
     & ISTUV_G,IENUV_G,IF_UV_G,ISTSC_G,IENSC_G,IF_SCALARS_G,&
     & ISTUV,IENUV,IF_UV,ISTSC,IENSC,IF_SCALARS)

    IF_FS = 2*IF_UV + IF_SCALARS
    IF_GP = 2*IF_UV_G+IF_SCALARS_G
    DO JFLD=1,IF_UV_G
      IPTRGP(JFLD) = ISHFUV_G(ISTUV_G+JFLD-1)
      IPTRGP(JFLD+IF_UV_G) = KF_UV_G+ISHFUV_G(ISTUV_G+JFLD-1)
    ENDDO
    DO JFLD=1,IF_SCALARS_G
      IPTRGP(JFLD+2*IF_UV_G) = 2*KF_UV_G+ISHFSC_G(ISTSC_G+JFLD-1)
    ENDDO
    DO JFLD=1,IF_UV
      IPTRSPUV(JFLD) = ISTUV+JFLD-1
    ENDDO
    DO JFLD=1,IF_SCALARS
      IPTRSPSC(JFLD) = ISTSC+JFLD-1
    ENDDO

    IF(IF_UV_G > 0 .AND. IF_SCALARS_G > 0) THEN
      CALL FTDIR_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,&
       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP,PGP=PGP)
    ELSEIF(IF_UV_G > 0) THEN
      CALL FTDIR_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,&
       & KVSETUV=IVSETUV(ISTUV_G:IENUV_G),&
       & KPTRGP=IPTRGP,PGP=PGP)
    ELSEIF(IF_SCALARS_G > 0) THEN
      CALL FTDIR_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,&
       & KVSETSC=IVSETSC(ISTSC_G:IENSC_G),KPTRGP=IPTRGP,PGP=PGP)
    ENDIF

    !$ACC DATA COPYIN(FOUBUF_IN)

    CALL LTDIR_CTL(IF_FS,IF_UV,IF_SCALARS, &
     & PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
     & KFLDPTRUV=IPTRSPUV,KFLDPTRSC=IPTRSPSC)

    !$ACC END DATA
  ENDDO
ELSE
  ! No splitting of fields, transform done in one go
  if (present(pgpuv)) then
    do jf=1,kf_uv_g
      write(nout,*) "pgpuv:",jf,mnx3(pgpuv(:,:,jf,:))
    end do
  end if

  if (present(pgp)) then
    do jf=1,kf_scalars_g
      write(nout,*) "pgp:",jf,mnx2(pgp(:,jf,:))
    end do
  end if

  if (present(pgp2)) then
    do jf=1,size(pgp2,2)
      write(nout,*) "pgp2:",jf,mnx2(pgp2(:,jf,:))
    end do
  end if

  if (present(pgp3a)) then
    do jf=1,size(pgp3a,3)
      write(nout,*) "pgp3a:",jf,mnx3(pgp3a(:,:,jf,:))
    end do
  end if

  if (present(pgp3b)) then
    do jf=1,size(pgp3b,3)
      write(nout,*) "pgp3b:",jf,mnx3(pgp3b(:,:,jf,:))
    end do
  end if

  CALL FTDIR_CTL(KF_UV_G,KF_SCALARS_G,KF_GP,KF_FS,KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
   & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
   & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2)

   CALL LTDIR_CTL(KF_FS,KF_UV,KF_SCALARS,PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
    & PSPSC3A=PSPSC3A,PSPSC3B=PSPSC3B,PSPSC2=PSPSC2)

  write(nout,*) "vor/div:",present(pspvor),present(pspdiv)
  if (kf_uv > 0) then
    do jf=1,kf_uv
      write(nout,*) "vor:",jf,minval(pspvor(jf,:)),maxval(pspvor(jf,:))
      write(nout,*) "div:",jf,minval(pspdiv(jf,:)),maxval(pspdiv(jf,:))
    end do
  end if

  if (kf_scalars) then
    if (present(pspscalar)) then
      do jf=1,kf_scalars
        write(nout,*) "scalar:",jf,minval(pspscalar(jf,:)),maxval(pspscalar(jf,:))
      end do
    else
      if (present(pspsc2).and.nf_sc2 > 0) then
      do jf=1,nf_sc2
        write(nout,*) "sc2:",jf,minval(pspsc2(jf,:)),maxval(pspsc2(jf,:))
      end do
      end if
      if (present(pspsc3a).and.nf_sc3a > 0) then
      do jf=1,nf_sc3a
        write(nout,*) "sc3a (=t):",jf,minval(pspsc3a(jf,:,1)),maxval(pspsc3a(jf,:,1))
      end do
      end if
    end if
  end if
ENDIF

END SUBROUTINE DIR_TRANS_CTL

function mnx2(x) result(mnx)
   real(kind=jprb),intent(in) :: x(:,:)

   real :: mnx(3)

   mnx(:) = (/minval(x),sum(x)/size(x),maxval(x)/)
end function

function mnx3(x) result(mnx)
   real(kind=jprb),intent(in) :: x(:,:,:)

   real :: mnx(3)

   mnx(:) = (/minval(x),sum(x)/size(x),maxval(x)/)
end function

END MODULE DIR_TRANS_CTL_MOD
