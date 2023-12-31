! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE DIR_TRANS(PSPVOR,PSPDIV,PSPSCALAR,PSPSC3A,PSPSC3B,PSPSC2,&
& LDLATLON,KPROMA,KVSETUV,KVSETSC,KRESOL,KVSETSC3A,KVSETSC3B,KVSETSC2,&
& PGP,PGPUV,PGP3A,PGP3B,PGP2)


!**** *DIR_TRANS* - Direct spectral transform (from grid-point to spectral).

!     Purpose.
!     --------
!        Interface routine for the direct spectral transform

!**   Interface.
!     ----------
!     CALL DIR_TRANS(...)

!     Explicit arguments : All arguments except from PGP are optional.
!     --------------------
!     PSPVOR(:,:) - spectral vorticity (output)
!     PSPDIV(:,:) - spectral divergence (output)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (output)
!     PSPSC3A(:,:,:) - alternative to use of PSPSCALAR, see PGP3A below (input)
!     PSPSC3B(:,:,:) - alternative to use of PSPSCALAR, see PGP3B below (input)
!     PSPSC2(:,:)  - alternative to use of PSPSCALAR, see PGP2 below (input)
!     LDLATLON   - indicating if regular lat-lon is the input data
!     KPROMA      - required blocking factor for gridpoint output
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
!     KVSETSC3A(:) - as KVESETSC for PSPSC3A (distribution on first dimension)
!     KVSETSC3B(:) - as KVESETSC for PSPSC3B (distribution on first dimension)
!     KVSETSC2(:) - as KVESETSC for PSPSC2 (distribution on first dimension)
!     KRESOL   - resolution tag  which is required ,default is the
!                first defined resulution (input)
!     PGP(:,:,:) - gridpoint fields (input)
!                  PGP need to  dimensioned (NPROMA,IF_GP,NGPBLKS) where
!                  NPROMA is the blocking factor, IF_GP the total number
!                  of output fields and NGPBLKS the number of NPROMA blocks.
!                  The ordering of the output fields is as follows (all
!                  parts are optional depending on the input switches):
!
!     u             : IF_UV_G fields (if psvor present)
!     v             : IF_UV_G fields (if psvor present)
!     scalar fields : IF_SCALARS_G fields (if pspscalar present)
!
!     Here IF_UV_G is the GLOBAL number of u/v fields as given by the length
!     of KVSETUV (or by PSPVOR if no split in spectral 'b-set' direction
!     IF_SCALARS_G is the GLOBAL number of scalar fields as giben by the
!     length of KVESETSC (or by number of fields in PSPSCALAR if no spectral
!     'b-set' split

!     As an alternative to using PGP you can also use a combination of the
!     following arrays. The reason for introducing these alternative ways
!     of calling DIR_TRANS is to avoid uneccessary copies where your data
!     structures don't fit in to the 'PSPVOR,PSPDIV, PSPSCALAR, PGP' layout.
!     The use of any of these precludes the use of PGP and vice versa.

!     PGPUV(:,:,:,:) - the 'u-v' related grid-point variables in the order
!                      described for PGP. The second dimension of PGPUV should
!                      be the same as the "global" first dimension of
!                      PSPVOR,PSPDIV (in the IFS this is the number of levels)
!                      PGPUV need to be dimensioned(NPROMA,ILEVS,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (u,v)
!     PGP3A(:,:,:,:) - grid-point array directly connected with PSPSC3A
!                      dimensioned(NPROMA,ILEVS,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (the same as in
!                      PSPSC3A )
!     PGP3B(:,:,:,:) - grid-point array directly connected with PSPSC3B
!                      dimensioned(NPROMA,ILEVS,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (the same as in
!                      PSPSC3B)
!     PGP2(:,:,:)    - grid-point array directly connected with PSPSC2
!                      dimensioned(NPROMA,IFLDS,NGPBLKS)
!                      IFLDS is the number of 'variables' (the same as in
!                      PSPSC2 )
!
!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  LTDIR_CTL   - control of Legendre transform
!                 FTDIR_CTL   - control of Fourier transform

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE TPM_GEN,ONLY : NERR, NOUT,NPROMATR
USE TPM_TRANS,ONLY : LDIVGP,LSCDERS,LUVDER,LVORGP,LATLON,NF_SC2,NF_SC3A,NF_SC3B,&
  NGPBLKS,NPROMA
USE TPM_DISTR       ,ONLY : D, NPRTRV, MYSETV,myproc
USE TPM_FIELDS      ,ONLY : IF_FS_DIR,IF_FS_DIR0,NFLEV,NFLEV0,DTDZBA,DTDZBS,DTDZCA,DTDZCS
USE TPM_FLT, ONLY: S
USE TPM_GEOMETRY ,ONLY : G
USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE DIR_TRANS_CTL_MOD ,ONLY : DIR_TRANS_CTL
USE LTDIR_CTL_MOD   ,ONLY : LTDIR_CTL
USE FTDIR_CTL_MOD   ,ONLY : FTDIR_CTL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK,  JPHOOK

IMPLICIT NONE

REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPVOR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPDIV(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSCALAR(:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3A(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC3B(:,:,:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PSPSC2(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KPROMA
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETUV(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3A(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC3B(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KVSETSC2(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KRESOL
LOGICAL   ,OPTIONAL, INTENT(IN) :: LDLATLON
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGP(:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGPUV(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGP3A(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGP3B(:,:,:,:)
REAL(KIND=JPRB),OPTIONAL    ,INTENT(IN) :: PGP2(:,:,:)

INTEGER(KIND=JPIM) :: IUBOUND(4),J
INTEGER(KIND=JPIM) :: IF_UV,IF_UV_G,IF_SCALARS,IF_SCALARS_G,IF_FS,IF_GP
INTEGER(KIND=JPIM) :: IF_SC2_G,IF_SC3A_G,IF_SC3B_G
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DIR_TRANS',0,ZHOOK_HANDLE)

CALL SET_RESOL(KRESOL)

! Set defaults
IF_UV = 0
IF_UV_G = 0
IF_SCALARS = 0
IF_SCALARS_G = 0
NF_SC2 = 0
NF_SC3A = 0
NF_SC3B = 0
IF_SC2_G = 0
IF_SC3A_G = 0
IF_SC3B_G = 0
NPROMA = D%NGPTOT
! This is for use in TRGTOL which is shared with adjoint inverse transform
LSCDERS=.FALSE.
LVORGP=.FALSE.
LDIVGP=.FALSE.
LUVDER=.FALSE.
LATLON=.FALSE.

! Decide requirements
! GP uv_g=nflevg=dim1(spvor), scal_g=nf=dim1(spscalar), sc2_g=dim1(spscalar)
! GP sc3a_g=dim1(sp3a) sc3b_g=dim1(sp3b) scal_g=sc2g+sc3ag+sc3bg
! SP uv=my vsetuv, scal=my vsetsc, sc2=my vsetsc2, sc3a/b=my vsetsc3a/b

IF(PRESENT(KVSETUV)) THEN
  IF_UV_G = UBOUND(KVSETUV,1)
  IF(ANY(KVSETUV(1:IF_UV_G) > NPRTRV .OR. KVSETUV(1:IF_UV_G) < 1)) THEN
    CALL ABORT_TRANS('DIR_TRANS:KVSETUV TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
  ENDIF

  IF_UV = COUNT(KVSETUV(1:IF_UV_G) == MYSETV)
ELSEIF(PRESENT(PSPVOR)) THEN
  IF_UV = UBOUND(PSPVOR,1)
  IF_UV_G = IF_UV
ENDIF

IF(PRESENT(KVSETSC)) THEN
  IF_SCALARS_G = UBOUND(KVSETSC,1)
  IF(ANY(KVSETSC(1:IF_SCALARS_G) > NPRTRV .OR. KVSETSC(1:IF_SCALARS_G) < 1)) THEN
    CALL ABORT_TRANS('DIR_TRANS:KVSETSC TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
  ENDIF

  IF_SCALARS = COUNT(KVSETSC(1:IF_SCALARS_G) == MYSETV)
ELSEIF(PRESENT(PSPSCALAR)) THEN
  IF_SCALARS = UBOUND(PSPSCALAR,1)
  IF_SCALARS_G = IF_SCALARS
ENDIF

IF(PRESENT(KVSETSC2)) THEN
  IF(.NOT.PRESENT(PSPSC2)) CALL ABORT_TRANS('DIR_TRANS:KVSETSC2 BUT NOT PSPSC2')

  IF_SC2_G = UBOUND(KVSETSC2,1)
  IF_SCALARS_G = IF_SCALARS_G+IF_SC2_G
  DO J=1,UBOUND(KVSETSC2,1)
    IF(KVSETSC2(J) > NPRTRV .OR. KVSETSC2(J) < 1) THEN
      CALL ABORT_TRANS('DIR_TRANS:KVSETSC2 TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSETSC2(J) == MYSETV) THEN
      IF_SCALARS = IF_SCALARS+1
      NF_SC2 = NF_SC2+1
    ENDIF
  ENDDO
ELSEIF(PRESENT(PSPSC2)) THEN
  IF_SC2_G = UBOUND(PSPSC2,1)
  NF_SC2   = UBOUND(PSPSC2,1)
  IF_SCALARS = IF_SCALARS+NF_SC2
  IF_SCALARS_G = IF_SCALARS_G +IF_SC2_G
ENDIF

IF(PRESENT(KVSETSC3A)) THEN
  IF(.NOT.PRESENT(PSPSC3A)) CALL ABORT_TRANS('DIR_TRANS:KVSETSC3A BUT NOT PSPSC3A')

  IF_SC3A_G = UBOUND(KVSETSC3A,1)
  IF_SCALARS_G = IF_SCALARS_G+IF_SC3A_G*UBOUND(PSPSC3A,3)
  DO J=1,UBOUND(KVSETSC3A,1)
    IF(KVSETSC3A(J) > NPRTRV .OR. KVSETSC3A(J) < 1) THEN
      CALL ABORT_TRANS('DIR_TRANS:KVSETSC3A TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSETSC3A(J) == MYSETV) THEN
      IF_SCALARS = IF_SCALARS+UBOUND(PSPSC3A,3)
      NF_SC3A = NF_SC3A+1
    ENDIF
  ENDDO
ELSEIF(PRESENT(PSPSC3A)) THEN
  IF_SCALARS = IF_SCALARS+UBOUND(PSPSC3A,1)*UBOUND(PSPSC3A,3)
  IF_SC3A_G = UBOUND(PSPSC3A,1)
  IF_SCALARS_G = IF_SCALARS_G +IF_SC3A_G*UBOUND(PSPSC3A,3)
  NF_SC3A = UBOUND(PSPSC3A,1)
ENDIF

IF(PRESENT(KVSETSC3B)) THEN
  IF(.NOT.PRESENT(PSPSC3B)) THEN
    CALL ABORT_TRANS('DIR_TRANS:KVSETSC3B BUT NOT PSPSC3B')
  ENDIF
  IF_SC3B_G = UBOUND(KVSETSC3B,1)
  IF_SCALARS_G = IF_SCALARS_G+IF_SC3B_G*UBOUND(PSPSC3B,3)
  DO J=1,UBOUND(KVSETSC3B,1)
    IF(KVSETSC3B(J) > NPRTRV .OR. KVSETSC3B(J) < 1) THEN
      CALL ABORT_TRANS('DIR_TRANS:KVSETSC3B TOO LONG OR CONTAINS VALUES OUTSIDE RANGE')
    ENDIF
    IF(KVSETSC3B(J) == MYSETV) THEN
      IF_SCALARS = IF_SCALARS+UBOUND(PSPSC3B,3)
      NF_SC3B = NF_SC3B+1
    ENDIF
  ENDDO
ELSEIF(PRESENT(PSPSC3B)) THEN
  IF_SCALARS = IF_SCALARS+UBOUND(PSPSC3B,1)*UBOUND(PSPSC3B,3)
  IF_SC3B_G = UBOUND(PSPSC3B,1)
  IF_SCALARS_G = IF_SCALARS_G +IF_SC3B_G*UBOUND(PSPSC3B,3)
  NF_SC3B = UBOUND(PSPSC3B,1)
ENDIF

IF(PRESENT(KPROMA)) NPROMA = KPROMA
IF(PRESENT(LDLATLON)) LATLON = LDLATLON

NGPBLKS = (D%NGPTOT-1)/NPROMA+1

IF_FS = 2*IF_UV+IF_SCALARS
!D%IADJUST_D=0
!IF(MOD(IF_FS,2)==1) THEN
!  IF_FS = IF_FS + 1
!  D%IADJUST_D=1
!ENDIF

IF_GP = 2*IF_UV_G+IF_SCALARS_G

write(nout,*) "UV fields:",IF_UV_G
write(nout,*) "SC fields:",IF_SCALARS_G
write(nout,*) "SC2 fields:",IF_SC2_G
write(nout,*) "SC3A fields:",IF_SC3A_G
write(nout,*) "SC3B fields:",IF_SC3B_G
write(nout,*) "GP fields (2*UV_G+SCALAR_G):",IF_GP
write(nout,*) "SP fields (2*UV+SCALAR):",IF_FS

! add additional post-processing requirements
! (copied from setup_trans.F90. Or does this need to be different here than in setup_trans.F90?)
!IF_PP = 2*NFLEV
!IF_PP = 0

! How do I get the current number of levels? For now I use: (Andreas)
!NFLEV = NFLEV0

! set currently used array sizes for the GPU arrays: 
IF_FS_DIR=2*IF_FS+2!2*(2*IF_UV+NFLEV+2+IF_PP)
print*,"dir_trans: IF_FS_DIR=",IF_FS_DIR," IF_FS_DIR0=",IF_FS_DIR0

DTDZBA=IF_FS_DIR
DTDZBS=IF_FS_DIR
DTDZCA=IF_FS_DIR
DTDZCS=IF_FS_DIR

IF (IF_UV > 0) THEN
  IF(.NOT. PRESENT(PSPVOR) ) THEN
    CALL ABORT_TRANS('DIR_TRANS : IF_UV > 0 BUT PSPVOR MISSING')
  ENDIF
  IF(UBOUND(PSPVOR,1) < IF_UV) THEN
    CALL ABORT_TRANS('DIR_TRANS : PSPVOR TOO SHORT')
  ENDIF
  IF(.NOT. PRESENT(PSPDIV) ) THEN
    CALL ABORT_TRANS('DIR_TRANS : PSPVOR PRESENT BUT PSPDIV MISSING')
  ENDIF
  IF(UBOUND(PSPDIV,1) /= IF_UV) THEN
    CALL ABORT_TRANS('DIR_TRANS : INCONSISTENT FIRST DIM. OF PSPVOR AND PSPDIV')
  ENDIF
ENDIF

IF (IF_SCALARS > 0) THEN
  IF(PRESENT(PSPSCALAR)) THEN
    IF(UBOUND(PSPSCALAR,1) < IF_SCALARS) THEN
      CALL ABORT_TRANS('DIR_TRANS : PSPSCALAR TOO SHORT')
    ENDIF
    IF(PRESENT(PSPSC3A))THEN
      CALL ABORT_TRANS('DIR_TRANS : PSPSCALAR AND PSPSC3A BOTH PRESENT')
    ENDIF
    IF(PRESENT(PSPSC3B))THEN
      CALL ABORT_TRANS('DIR_TRANS : PSPSCALAR AND PSPSC3B BOTH PRESENT')
    ENDIF
    IF(PRESENT(PSPSC2))THEN
      CALL ABORT_TRANS('DIR_TRANS : PSPSCALAR AND PSPSC2 BOTH PRESENT')
    ENDIF
  ENDIF
ENDIF

IF(NPRTRV >1) THEN
  IF(IF_UV > 0 .AND. .NOT. PRESENT(KVSETUV)) THEN
    CALL ABORT_TRANS('DIR_TRANS: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
  IF(PRESENT(PSPSCALAR) .AND. .NOT. PRESENT(KVSETSC)) THEN
    CALL ABORT_TRANS('DIR_TRANS: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
  IF(PRESENT(PSPSC2) .AND. .NOT. PRESENT(KVSETSC2)) THEN
    CALL ABORT_TRANS('DIR_TRANS: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
  IF(PRESENT(PSPSC3A) .AND. .NOT. PRESENT(KVSETSC3A)) THEN
    CALL ABORT_TRANS('DIR_TRANS: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
  IF(PRESENT(PSPSC3B) .AND. .NOT. PRESENT(KVSETSC3B)) THEN
    CALL ABORT_TRANS('DIR_TRANS: SPECIFY VERTICAL SPECTRAL DISTRIBUTION!')
  ENDIF
ENDIF

IF(PRESENT(PGP)) THEN
  IUBOUND(1:3)=UBOUND(PGP)
  IF(IUBOUND(1) < NPROMA) THEN
    CALL ABORT_TRANS('DIR_TRANS:FIRST DIMENSION OF PGP TOO SMALL ')
  ENDIF
  IF(IUBOUND(2) < IF_GP) THEN
    CALL ABORT_TRANS('DIR_TRANS:SECOND DIMENSION OF PGP TOO SMALL ')
  ENDIF
  IF(IUBOUND(3) < NGPBLKS) THEN
    CALL ABORT_TRANS('DIR_TRANS:THIRD DIMENSION OF PGP TOO SMALL ')
  ENDIF
ENDIF

IF(PRESENT(PGPUV)) THEN
  IF(.NOT.PRESENT(PSPVOR)) THEN
    CALL ABORT_TRANS('DIR_TRANS:PSPVOR HAS TO BE PRESENT WHEN PGPUV IS')
  ENDIF
  IUBOUND=UBOUND(PGPUV)
  IF(IUBOUND(1) < NPROMA) THEN
    CALL ABORT_TRANS('DIR_TRANS:FIRST DIMENSION OF PGPUV TOO SMALL ')
  ENDIF
  IF(IUBOUND(2) /= IF_UV_G) THEN
    CALL ABORT_TRANS('DIR_TRANS:SEC. DIMENSION OF PGPUV INCONSISTENT ')
  ENDIF
  IF(IUBOUND(3) < 2) THEN
    CALL ABORT_TRANS('DIR_TRANS:THIRD DIMENSION OF PGPUV TOO SMALL ')
  ENDIF
  IF(IUBOUND(4) < NGPBLKS) THEN
    CALL ABORT_TRANS('DIR_TRANS:FOURTH DIMENSION OF PGPUV TOO SMALL ')
  ENDIF
ENDIF

IF(PRESENT(PGP2)) THEN
  IF(.NOT.PRESENT(PSPSC2)) THEN
    CALL ABORT_TRANS('DIR_TRANS:PSPSC2 HAS TO BE PRESENT WHEN PGP2 IS')
  ENDIF
ENDIF
IF(IF_SC2_G > 0) THEN
  IF(PRESENT(PGP2)) THEN
    IUBOUND(1:3)=UBOUND(PGP2)
    IF(IUBOUND(1) < NPROMA) THEN
      CALL ABORT_TRANS('DIR_TRANS:FIRST DIMENSION OF PGP2 TOO SMALL ')
    ENDIF
    IF(IUBOUND(2) /= IF_SC2_G) THEN
      CALL ABORT_TRANS('DIR_TRANS:SEC. DIMENSION OF PGP2 INCONSISTENT')
    ENDIF
    IF(IUBOUND(3) < NGPBLKS) THEN
      CALL ABORT_TRANS('DIR_TRANS:THIRD DIMENSION OF PGP2 TOO SMALL ')
    ENDIF
  ELSE
    CALL ABORT_TRANS('DIR_TRANS:PGP2 MISSING')
  ENDIF
ENDIF

IF(PRESENT(PGP3A)) THEN
  IF(.NOT.PRESENT(PSPSC3A)) THEN
    CALL ABORT_TRANS('DIR_TRANS:PSPSC3A HAS TO BE PRESENT WHEN PGP3A IS')
  ENDIF
ENDIF
IF(IF_SC3A_G > 0) THEN
  IF(PRESENT(PGP3A)) THEN
    IUBOUND=UBOUND(PGP3A)
    IF(IUBOUND(1) < NPROMA) THEN
      CALL ABORT_TRANS('DIR_TRANS:FIRST DIMENSION OF PGP3A TOO SMALL ')
    ENDIF
    IF(IUBOUND(2) /= IF_SC3A_G) THEN
      CALL ABORT_TRANS('DIR_TRANS:SEC. DIMENSION OF PGP3A INCONSISTENT ')
    ENDIF
    IF(IUBOUND(3) /= UBOUND(PSPSC3A,3) ) THEN
      CALL ABORT_TRANS('DIR_TRANS:THIRD DIMENSION OF PGP3A INCONSISTENT ')
    ENDIF
    IF(IUBOUND(4) < NGPBLKS) THEN
      CALL ABORT_TRANS('DIR_TRANS:FOURTH DIMENSION OF PGP3A TOO SMALL ')
    ENDIF
  ELSE
    CALL ABORT_TRANS('DIR_TRANS:PGP3A MISSING')
  ENDIF
ENDIF

IF(PRESENT(PGP3B)) THEN
  IF(.NOT.PRESENT(PSPSC3B)) THEN
    CALL ABORT_TRANS('DIR_TRANS:PSPSC3B HAS TO BE PRESENT WHEN PGP3B IS')
  ENDIF
ENDIF
IF(IF_SC3B_G > 0) THEN
  IF(PRESENT(PGP3B)) THEN
    IUBOUND=UBOUND(PGP3B)
    IF(IUBOUND(1) < NPROMA) THEN
      CALL ABORT_TRANS('DIR_TRANS:FIRST DIMENSION OF PGP3B TOO SMALL ')
    ENDIF
    IF(IUBOUND(2) /= IF_SC3B_G) THEN
      CALL ABORT_TRANS('DIR_TRANS:SEC. DIMENSION OF PGP3B INCONSISTENT ')
    ENDIF
    IF(IUBOUND(3) /= UBOUND(PSPSC3B,3) ) THEN
      CALL ABORT_TRANS('DIR_TRANS:THIRD DIMENSION OF PGP3B INCONSISTENT ')
    ENDIF
    IF(IUBOUND(4) < NGPBLKS) THEN
      CALL ABORT_TRANS('DIR_TRANS:FOURTH DIMENSION OF PGP3B TOO SMALL ')
    ENDIF
  ELSE
    CALL ABORT_TRANS('DIR_TRANS:PGP3B MISSING')
  ENDIF
ENDIF

IF (NPROMATR > 0 .AND. IF_GP > NPROMATR) THEN
  CALL DIR_TRANS_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,IF_UV,IF_SCALARS,&
    PSPVOR,PSPDIV,PSPSCALAR,KVSETUV,KVSETSC,PGP,&
    PSPSC3A,PSPSC3B,PSPSC2,KVSETSC3A,KVSETSC3B,KVSETSC2,PGPUV,PGP3A,PGP3B,PGP2)
ELSE
  CALL FTDIR_CTL(IF_UV_G,IF_SCALARS_G,IF_GP,IF_FS,KVSETUV=KVSETUV,KVSETSC=KVSETSC,&
   & KVSETSC3A=KVSETSC3A,KVSETSC3B=KVSETSC3B,KVSETSC2=KVSETSC2,&
   & PGP=PGP,PGPUV=PGPUV,PGP3A=PGP3A,PGP3B=PGP3B,PGP2=PGP2)

  CALL LTDIR_CTL(IF_FS,IF_UV,IF_SCALARS,PSPVOR=PSPVOR,PSPDIV=PSPDIV,PSPSCALAR=PSPSCALAR,&
   & PSPSC3A=PSPSC3A,PSPSC3B=PSPSC3B,PSPSC2=PSPSC2)
ENDIF

IF (LHOOK) CALL DR_HOOK('DIR_TRANS',1,ZHOOK_HANDLE)
END SUBROUTINE DIR_TRANS

