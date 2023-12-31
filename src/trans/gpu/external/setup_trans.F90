! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SETUP_TRANS(KSMAX,KDGL,KDLON,KLOEN,LDSPLIT,PSTRET,&
&KFLEV,KTMAX,KRESOL,PWEIGHT,LDGRIDONLY,LDUSERPNM,LDKEEPRPNM,LDUSEFLT,&
&LDSPSETUPONLY,LDPNMONLY,LDUSEFFTW,&
&LDLL,LDSHIFTLL,CDIO_LEGPOL,CDLEGPOLFNAME,KLEGPOLPTR,KLEGPOLPTR_LEN)

!**** *SETUP_TRANS* - Setup transform package for specific resolution

!     Purpose.
!     --------
!     To setup for making spectral transforms. Each call to this routine
!     creates a new resolution up to a maximum of NMAX_RESOL set up in
!     SETUP_TRANS0. You need to call SETUP_TRANS0 before this routine can
!     be called.

!**   Interface.
!     ----------
!     CALL SETUP_TRANS(...)

!     Explicit arguments : KLOEN,LDSPLIT are optional arguments
!     --------------------
!     KSMAX - spectral truncation required
!     KDGL  - number of Gaussian latitudes
!     KDLON - number of points on each Gaussian latitude [2*KDGL]
!     KLOEN(:) - number of points on each Gaussian latitude [2*KDGL]
!     LDSPLIT - true if split latitudes in grid-point space [false]
!     KTMAX - truncation order for tendencies?
!     KRESOL - the resolution identifier
!     PWEIGHT - the weight per grid-point (for a weighted distribution)
!     LDGRIDONLY - true if only grid space is required

!     KSMAX,KDGL,KTMAX and KLOEN are GLOBAL variables desribing the resolution
!     in spectral and grid-point space

!     LDSPLIT describe the distribution among processors of grid-point data and
!     has no relevance if you are using a single processor

!     PSTRET     - stretching factor - for the case the Legendre polynomials are
!                  computed on the stretched sphere - works with LSOUTHPNM
!     LDUSEFLT   - use Fast Legandre Transform (Butterfly algorithm)
!     LDUSERPNM  - Use Belusov algorithm to compute legendre pol. (else new alg.)
!     LDKEEPRPNM - Keep Legendre Polynomials (only applicable when using
!                  FLT, otherwise always kept)
!     LDPNMONLY  - Compute the Legendre polynomials only, not the FFTs.
!     LDUSEFFTW    - Use FFTW for FFTs
!     LDLL                 - Setup second set of input/output latitudes
!                                 the number of input/output latitudes to transform is equal KDGL
!                                 or KDGL+2 in the case that includes poles + equator
!                                 the number of input/output longitudes to transform is 2*KDGL
!     LDSHIFTLL       - Shift output lon/lat data by 0.5*dx and 0.5*dy
!     CDIO_LEGPOL  - IO option on Legendre polinomials :  N.B. Only works for NPROC=1
!                    Options:
!                    'READF' -  read Leg.Pol. from file CDLEGPOLFNAME
!                    'WRITEF' - write Leg.Pol. to file CDLEGPOLFNAME
!                    'MEMBUF' - Leg. Pol provided in shared memory segment pointed to by KLEGPOLPTR of
!                               length KLEGPOLPTR_LEN
!     CDLEGPOLFNAME - file name for Leg.Pol. IO
!     KLEGPOLPTR    - pointer to Legendre polynomials memory segment
!     KLEGPOLPTR_LEN  - length of  Legendre polynomials memory segment

!     Method.
!     -------

!     Externals.  SET_RESOL   - set resolution
!     ----------  SETUP_DIMS  - setup distribution independent dimensions
!                 SUMP_TRANS_PRELEG - first part of setup of distr. environment
!                 SULEG - Compute Legandre polonomial and Gaussian
!                         Latitudes and Weights
!                 SUMP_TRANS - Second part of setup of distributed environment
!                 SUFFT - setup for FFT
!                 SHAREDMEM_CREATE - create memory buffer for Leg.pol.

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        Daan Degrauwe : Mar 2012 E'-zone dimensions
!        R. El Khatib 09-Aug-2012 %LAM in GEOM_TYPE
!        R. El Khatib 14-Jun-2013 PSTRET, LDPNMONLY, LENABLED
!        G. Mozdzynski : Oct 2014 Support f
!        N. Wedi       : Apr 2015 Support dual set of lat/lon
!        G. Mozdzynski : Jun 2015 Support alternative FFTs to FFTW
!        M.Hamrud/W.Deconinck : July 2015 IO options for Legenndre polynomials
!        R. El Khatib 07-Mar-2016 Better flexibility for Legendre polynomials computation in stretched mode
!     ------------------------------------------------------------------

USE PARKIND1        ,ONLY : JPIM     ,JPRB ,  JPRD
USE PARKIND_ECTRANS ,ONLY : JPRBT
USE, INTRINSIC :: ISO_C_BINDING, ONLY:  C_PTR, C_INT,C_ASSOCIATED,C_SIZE_T

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NOUT, MSETUP0, NCUR_RESOL, NDEF_RESOL, &
     &                      NMAX_RESOL, NPRINTLEV, LENABLED, NERR
USE TPM_DIM         ,ONLY : R, DIM_RESOL, R_NSMAX,R_NTMAX, R_NDGNH, R_NDGL, R_NNOEXTZL
USE TPM_DISTR       ,ONLY : D, DISTR_RESOL,NPROC,nprtrv, D_NUMP,D_MYMS,D_NSTAGT0B,D_NSTAGT1B,D_NPROCL,D_NPNTGTB1, D_NASM0, &
& D_NSTAGTF,D_MSTABF,D_NPNTGTB0,D_NPROCM,D_NPTRLS,mysetv,mysetw, MYPROC
USE TPM_GEOMETRY    ,ONLY : G, GEOM_RESOL, G_NDGLU, G_NMEN, G_NMEN_MAX,G_NLOEN, G_NLOEN_MAX
USE TPM_FIELDS      ,ONLY : FIELDS_RESOL, F,F_RW, ZIA,ZEPSNM,ZSOA1,ZAOA1,ISTAN,ISTAS,ZAIA,ZOA1,ZOA2, &
& ZAA,ZAS,LDZAA,LDZAS,TDZAA,TDZAS,&
& IZBS,ILDZBA,ILDZBS,ITDZBA0,ITDZBS0,&
& IZCA,IZCS,IZCST,ILDZCA,ILDZCS,ITDZCA0,ITDZCS0,&
& DZBAT,DZBST,DLDZBA,DLDZBS,DTDZBA0,DTDZBS0,&
& DZCAT,DZCST,DLDZCA,DLDZCS,DTDZCA0,DTDZCS0,&
& IF_FS_INV0,IF_FS_DIR0,NFLEV0,ZAA0,DZBST0,DZCAT0,&
& ZAS0,DZCST0,KMLOC0
! IZBA,IZCAT
USE TPM_FFT         ,ONLY : T, FFT_RESOL
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, FFTW_RESOL
#endif
USE TPM_FFTC        ,ONLY : TC, FFTC_RESOL
USE TPM_FLT
USE TPM_TRANS       ,ONLY : FOUBUF_IN, FOUBUF, ZGTF, ZAVE, ZMINGL, ZMAXGL, ZMINGPN, ZMAXGPN
USE TPM_CTL

USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE SETUP_DIMS_MOD  ,ONLY : SETUP_DIMS
USE SUMP_TRANS_MOD  ,ONLY : SUMP_TRANS
USE SUMP_TRANS_PRELEG_MOD ,ONLY : SUMP_TRANS_PRELEG
USE SULEG_MOD       ,ONLY : SULEG
USE PRE_SULEG_MOD   ,ONLY : PRE_SULEG
USE SUFFT_MOD       ,ONLY : SUFFT
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE SHAREDMEM_MOD    ,ONLY : SHAREDMEM_CREATE
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK,  JPHOOK
USE CUDA_DEVICE_MOD
USE PREPSNM_MOD     ,ONLY : PREPSNM
#ifdef _OPENACC
use openacc
#endif

!endif INTERFACE

IMPLICIT NONE

! Dummy arguments

INTEGER(KIND=JPIM) ,INTENT(IN) :: KSMAX,KDGL
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KDLON
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KLOEN(:)
LOGICAL            ,OPTIONAL,INTENT(IN) :: LDSPLIT
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KTMAX
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT):: KRESOL
REAL(KIND=JPRB)    ,OPTIONAL,INTENT(IN) :: PWEIGHT(:)
REAL(KIND=JPRB)    ,OPTIONAL,INTENT(IN) :: PSTRET
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN) :: KFLEV
LOGICAL   ,OPTIONAL,INTENT(IN):: LDGRIDONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFLT
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSERPNM
LOGICAL   ,OPTIONAL,INTENT(IN):: LDKEEPRPNM
LOGICAL   ,OPTIONAL,INTENT(IN):: LDSPSETUPONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDPNMONLY
LOGICAL   ,OPTIONAL,INTENT(IN):: LDUSEFFTW
LOGICAL   ,OPTIONAL,INTENT(IN):: LDLL
LOGICAL   ,OPTIONAL,INTENT(IN):: LDSHIFTLL
CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: CDIO_LEGPOL
CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: CDLEGPOLFNAME
TYPE(C_PTR) ,OPTIONAL,INTENT(IN) :: KLEGPOLPTR
INTEGER(C_SIZE_T) ,OPTIONAL,INTENT(IN) :: KLEGPOLPTR_LEN

!ifndef INTERFACE

! Local variables
INTEGER(KIND=JPIM),PARAMETER :: IMAXFLD=240
INTEGER(KIND=JPIM) :: JGL,JRES,IDEF_RESOL
INTEGER(KIND=JPIM) :: NFLEVL, JMLOC, KM, ILA, ILS, KMLOC, KDGLU, JK, i, J, IF_FS, IF_OUT_LT, IF_UV, IF_SCALARS
INTEGER(KIND=JPIM) :: IPPNUM, IF_PP, IF_FOUBUF

LOGICAL :: LLP1,LLP2, LLSPSETUPONLY
REAL(KIND=JPRD)    :: ZTIME0,ZTIME1,ZTIME2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

integer :: idevtype, inumdevs, mygpu, iunit, istat, idev

#include "user_clock.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SETUP_TRANS',0,ZHOOK_HANDLE)

IF (MSETUP0 == 0) CALL ABORT_TRANS('SETUP_TRANS0 HAS NOT BEEN CALLED')

LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1
IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SETUP_TRANS ==='

! Allocate resolution dependent structures
IF(.NOT. ALLOCATED(DIM_RESOL)) THEN
  IDEF_RESOL = 1
  ALLOCATE(DIM_RESOL(NMAX_RESOL))
  ALLOCATE(FIELDS_RESOL(NMAX_RESOL))
  ALLOCATE(GEOM_RESOL(NMAX_RESOL))
  ALLOCATE(DISTR_RESOL(NMAX_RESOL))
  ALLOCATE(FFT_RESOL(NMAX_RESOL))
#ifdef WITH_FFTW
  ALLOCATE(FFTW_RESOL(NMAX_RESOL))
#endif
  ALLOCATE(FFTC_RESOL(NMAX_RESOL))
  ALLOCATE(FLT_RESOL(NMAX_RESOL))
  ALLOCATE(CTL_RESOL(NMAX_RESOL))
  GEOM_RESOL(:)%LAM=.FALSE.
  ALLOCATE(LENABLED(NMAX_RESOL))
  LENABLED(:)=.FALSE.
ELSE
  IDEF_RESOL = NMAX_RESOL+1
  DO JRES=1,NMAX_RESOL
    IF(.NOT.LENABLED(JRES)) THEN
      IDEF_RESOL = JRES
      EXIT
    ENDIF
  ENDDO
  IF(IDEF_RESOL > NMAX_RESOL) THEN
    CALL ABORT_TRANS('SETUP_TRANS:IDEF_RESOL > NMAX_RESOL')
  ENDIF
ENDIF

IF (PRESENT(KRESOL)) THEN
  KRESOL=IDEF_RESOL
ENDIF

! Point at structures due to be initialized
CALL SET_RESOL(IDEF_RESOL,LDSETUP=.TRUE.)

IF(LLP1) WRITE(NOUT,*) '=== DEFINING RESOLUTION ',NCUR_RESOL



! Defaults for optional arguments


G%LREDUCED_GRID = .FALSE.
G%RSTRET=1.0_JPRBT
D%LGRIDONLY = .FALSE.
D%LSPLIT = .FALSE.
D%LCPNMONLY=.FALSE.
S%LUSE_BELUSOV=.TRUE. ! use Belusov algorithm to compute RPNM array instead of per m
S%LKEEPRPNM=.FALSE. ! Keep Legendre polonomials (RPNM)
S%LUSEFLT=.FALSE. ! Use fast legendre transforms
#ifdef WITH_FFTW
TW%LFFTW=.FALSE. ! Use FFTW interface for FFTs
#endif
LLSPSETUPONLY = .FALSE. ! Only create distributed spectral setup
S%LDLL = .FALSE. ! use mapping to/from second set of latitudes
S%LSHIFTLL = .FALSE. ! shift output lat-lon by 0.5dx, 0.5dy
C%LREAD_LEGPOL = .FALSE.
C%LWRITE_LEGPOL = .FALSE.


! NON-OPTIONAL ARGUMENTS
R%NSMAX = KSMAX
R%NDGL  = KDGL
! E'-defaults
R%NNOEXTZL=0
R%NNOEXTZG=0

! IMPLICIT argument :
G%LAM = .FALSE.

IF(PRESENT(KDLON)) THEN
  R%NDLON = KDLON
ELSE
  R%NDLON = 2*R%NDGL
ENDIF

IF(PRESENT(LDLL)) THEN
  S%LDLL=LDLL
  IF( LDLL ) THEN
    S%NDLON=R%NDLON
    ! account for pole + equator
    R%NDGL=R%NDGL+2
    IF(PRESENT(LDSHIFTLL)) THEN
      S%LSHIFTLL = LDSHIFTLL
      ! geophysical (shifted) lat-lon without pole and equator
      IF(S%LSHIFTLL) R%NDGL=R%NDGL-2
    ENDIF
    S%NDGL=R%NDGL
  ENDIF
ENDIF

IF (R%NDGL <= 0 .OR. MOD(R%NDGL,2) /= 0) THEN
  CALL ABORT_TRANS ('SETUP_TRANS: KDGL IS NOT A POSITIVE, EVEN NUMBER')
ENDIF

! Optional arguments

ALLOCATE(G%NLOEN(R%NDGL))
IF(LLP2)WRITE(NOUT,9) 'NLOEN   ',SIZE(G%NLOEN   ),SHAPE(G%NLOEN   )
IF(PRESENT(KLOEN)) THEN
  IF( MINVAL(KLOEN(:)) <= 0 )THEN
     CALL ABORT_TRANS ('SETUP_TRANS: KLOEN INVALID (ONE or MORE POINTS <= 0)')
  ENDIF
  R%NDLON=MAXVAL(KLOEN(:))
  DO JGL=1,R%NDGL
    IF(KLOEN(JGL) /= R%NDLON) THEN
      G%LREDUCED_GRID = .TRUE.
      EXIT
    ENDIF
  ENDDO
ENDIF

IF (G%LREDUCED_GRID) THEN
  G%NLOEN(:) = KLOEN(1:R%NDGL)
ELSE
  G%NLOEN(:) = R%NDLON
ENDIF

IF(PRESENT(LDSPLIT)) THEN
  D%LSPLIT = LDSPLIT
ENDIF

IF(PRESENT(KTMAX)) THEN
  R%NTMAX = KTMAX
ELSE
  R%NTMAX = R%NSMAX
ENDIF

IF(PRESENT(PWEIGHT)) THEN
  D%LWEIGHTED_DISTR = .TRUE.
  IF( D%LWEIGHTED_DISTR .AND. .NOT.D%LSPLIT )THEN
    CALL ABORT_TRANS('SETUP_TRANS: LWEIGHTED_DISTR=T AND LSPLIT=F NOT SUPPORTED')
  ENDIF
  IF(SIZE(PWEIGHT) /= SUM(G%NLOEN(:)) )THEN
    CALL ABORT_TRANS('SETUP_TRANS:SIZE(PWEIGHT) /= SUM(G%NLOEN(:))')
  ENDIF
  IF( MINVAL(PWEIGHT(:)) < 0.0_JPRBT )THEN
    CALL ABORT_TRANS('SETUP_TRANS: INVALID WEIGHTS')
  ENDIF
  ALLOCATE(D%RWEIGHT(SIZE(PWEIGHT)))
  D%RWEIGHT(:)=PWEIGHT(:)
ELSE
  D%LWEIGHTED_DISTR = .FALSE.
ENDIF

IF(PRESENT(LDGRIDONLY)) THEN
  D%LGRIDONLY=LDGRIDONLY
ENDIF

IF(PRESENT(LDSPSETUPONLY)) THEN
  LLSPSETUPONLY=LDSPSETUPONLY
ENDIF

IF(PRESENT(LDPNMONLY)) THEN
  D%LCPNMONLY=LDPNMONLY
ENDIF


#ifdef WITH_FFTW
IF(PRESENT(LDUSEFFTW)) THEN
  TW%LFFTW=LDUSEFFTW
ENDIF
IF( LLSPSETUPONLY .OR. D%LGRIDONLY ) THEN
  TW%LFFTW = .FALSE.
ENDIF
#endif

S%LSOUTHPNM=.FALSE.
IF(PRESENT(PSTRET)) THEN
  IF (ABS(PSTRET-1.0_JPRBT)>100._JPRBT*EPSILON(1._JPRBT)) THEN
    G%RSTRET=PSTRET
    S%LSOUTHPNM=.TRUE.
  ENDIF
ENDIF

IF(PRESENT(CDIO_LEGPOL)) THEN
  IF(NPROC > 1) CALL  ABORT_TRANS('SETUP_TRANS:CDIO_LEGPOL OPTIONS ONLY FOR NPROC=1 ')
  IF(R%NSMAX > 511 ) S%LUSEFLT = .TRUE. !To save IO and memory
  IF(TRIM(CDIO_LEGPOL) == 'readf' .OR. TRIM(CDIO_LEGPOL) == 'READF' ) THEN
    IF(.NOT.PRESENT(CDLEGPOLFNAME)) CALL  ABORT_TRANS('SETUP_TRANS: CDLEGPOLFNAME ARGUMENT MISSING')
    C%LREAD_LEGPOL = .TRUE.
    C%CLEGPOLFNAME = TRIM(CDLEGPOLFNAME)
    C%CIO_TYPE='file'
  ELSEIF(TRIM(CDIO_LEGPOL) == 'writef' .OR. TRIM(CDIO_LEGPOL) == 'WRITEF') THEN
    IF(.NOT.PRESENT(CDLEGPOLFNAME)) CALL  ABORT_TRANS('SETUP_TRANS: CDLEGPOLFNAME ARGUMENT MISSING')
    C%LWRITE_LEGPOL = .TRUE.
    C%CLEGPOLFNAME = TRIM(CDLEGPOLFNAME)
    C%CIO_TYPE='file'
  ELSEIF(TRIM(CDIO_LEGPOL) == 'membuf' .OR. TRIM(CDIO_LEGPOL) == 'MEMBUF') THEN
    IF(.NOT.PRESENT(KLEGPOLPTR)) CALL  ABORT_TRANS('SETUP_TRANS: KLEGPOLPTR  ARGUMENT MISSING')
    IF(.NOT.C_ASSOCIATED(KLEGPOLPTR))  CALL  ABORT_TRANS('SETUP_TRANS: KLEGPOLPTR NULL POINTER')
    IF(.NOT.PRESENT(KLEGPOLPTR_LEN)) CALL  ABORT_TRANS('SETUP_TRANS: KLEGPOLPTR_LEN ARGUMENT MISSING')
    C%LREAD_LEGPOL = .TRUE.
    C%CIO_TYPE='mbuf'
    CALL SHAREDMEM_CREATE( C%STORAGE,KLEGPOLPTR,KLEGPOLPTR_LEN)
  ELSE
    WRITE(NERR,*) 'CDIO_LEGPOL ', TRIM(CDIO_LEGPOL)
    CALL  ABORT_TRANS('SETUP_TRANS:CDIO_LEGPOL UNKNOWN METHOD ')
  ENDIF
ENDIF

IF(PRESENT(LDUSEFLT)) THEN
  S%LUSEFLT=LDUSEFLT
ENDIF
IF(PRESENT(LDUSERPNM)) THEN
  S%LUSE_BELUSOV=LDUSERPNM
ENDIF
IF(PRESENT(LDKEEPRPNM)) THEN
  IF(S%LUSEFLT) THEN
    IF(LDKEEPRPNM.AND..NOT.LDUSERPNM) THEN
      CALL ABORT_TRANS('SETUP_TRANS: LDKEEPRPNM=true with LDUSERPNM=false')
    ENDIF
  ENDIF
  S%LKEEPRPNM=LDKEEPRPNM
ENDIF
!     Setup resolution dependent structures
!     -------------------------------------

write(nout,*) "Setup distribution independent dimensions"
CALL SETUP_DIMS

write(nout,*) "First part of setup of distributed environment"
CALL SUMP_TRANS_PRELEG

IF( .NOT.LLSPSETUPONLY ) THEN

  write(nout,*) "Compute Legendre polonomial and Gaussian Latitudes and Weights"
  CALL SULEG

  write(nout,*) "Second part of setup of distributed environment"
  CALL SUMP_TRANS
  CALL GSTATS(1802,0)

! Initialize Fast Fourier Transform package
  IF (.NOT.D%LCPNMONLY) CALL SUFFT
  CALL GSTATS(1802,1)
ELSE
  CALL PRE_SULEG
ENDIF

! Signal the current resolution is active
LENABLED(IDEF_RESOL)=.TRUE.
NDEF_RESOL = COUNT(LENABLED)

IF (LHOOK) CALL DR_HOOK('SETUP_TRANS',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

IF( .NOT.D%LGRIDONLY ) THEN

!allocating arrays for the GPU:
IF(PRESENT(KFLEV)) THEN
  NFLEV0 = KFLEV
!  NFLEVL = NFLEV0/NPRTRV
ELSE
  NFLEV0 = ceiling(REAL(IMAXFLD)/NPRTRV)
ENDIF

! need to get local rank to be able to set device (1GPU == 1 MPI-rank)
!ilocal_rank = 0
!call GETENV("OMPI_COMM_WORLD_LOCAL_RANK",comm_local_rank)
!read(comm_local_rank,'(I2)') ilocal_rank

iunit=300+myproc

! daand: removed this entirely
!#ifdef _OPENACC
!!!idevtype=acc_device_nvidia
!idevtype=acc_get_device_type()
!inumdevs = acc_get_num_devices(idevtype)
!mygpu = mod(MYPROC-1,inumdevs)
!CALL acc_set_device_num(mygpu, idevtype)
!mygpu = acc_get_device_num(idevtype)
!istat  = cuda_GetDevice(idev)
!write (0,*) 'cuda_getdevice returned code ',istat
!WRITE(0,*) '===now going to allocate GPU arrays on processor: ', myproc, ' device = ', mygpu, ' ',idev, ' of ', inumdevs
!#endif

!dimensions of matrices for Legendre Transforms for RAPS ?
!IF_OUT_LT = 5*NFLEV0+2
!IF_FS = 6*NFLEV0+3

! add additional post-processing requirements
!IF_PP = 2*NFLEV0
IF_PP = 0

! u/v + scalars 3d + scalars 2d
IF_UV = NFLEV0
! SCALARS INCLUDING DERIVATIVES
IF_SCALARS = NFLEV0 + 2*NFLEV0 + 1 + 2 + IF_PP
IF_OUT_LT = 4*IF_UV+3*NFLEV0+3+IF_PP
!IF_OUT_LT = 4*IF_UV+3*NFLEV0+3
!8*KF_UV+2*KF_SCALARS
!ILEI2 = 8*KF_UV + 2*KF_SCALARS + 2*KF_SCDERS
IF_FS_INV0=8*IF_UV+2*IF_SCALARS

! fields in Fourier space for inv trans the same
!IF_FS=4*IF_UV+1*NFLEV0+2
IF_FS=4*IF_UV+1*NFLEV0+2
! for derivatives u/v add
!IF_FS=IFS_FS+2*(2*NFLEV0)
! for each 3d scalar derivative add
IF_FS=IF_FS+2*NFLEV0 ! temperature
! for each 2d scalar derivative add
IF_FS=IF_FS+2 ! sfc pressure
IF_FS=IF_FS+IF_PP

! u/v + scalars for direct transforms
! plus postprocessing buffer
!ippnum=NFLEV0
IF_FS_DIR0=2*(2*IF_UV+NFLEV0+2+IF_PP)
!QUESTION: Why do we have NFLEV0 here? (Andreas)

! fields in Fourier space for dir trans
!IF_FS = 2*IF_UV + IF_SCALARS
! plus add 2*scalar_derivatives + add vorg/divg + 2*IF_UV for u/v zonal derivatives

write(nout,*)'setup_trans: if_uv=',if_uv,' if_out_lt=',if_out_lt,' IF_FS_DIR0=',IF_FS_DIR0,'IF_FS_INV0= ',IF_FS_INV0
IF(MOD(IF_FS,2)==1) IF_FS = IF_FS + 1

!leading and trailing dimensions of A for symmetric and antisymmetric cases
! (same for ltinv and ltdir)
LDZAA=R%NDGNH
LDZAS=R%NDGNH
TDZAA=(R%NTMAX+2)/2
TDZAS=(R%NTMAX+3)/2
print*,'R%NTMAX=',R%NTMAX
print*,'R%NSMAX=',R%NSMAX
!similarly for B (ltinv)
ILDZBA=(R%NSMAX+2)/2
ILDZBS=(R%NSMAX+3)/2
ITDZBA0=IF_FS_INV0
ITDZBS0=IF_FS_INV0

!similarly for C (ltinv)
ILDZCA=R%NDGNH
ILDZCS=R%NDGNH
ITDZCA0=IF_FS_INV0
ITDZCS0=IF_FS_INV0

!similarly for B (ltdir)
DLDZBA=R%NDGNH
DLDZBS=R%NDGNH
DTDZBA0=IF_FS_DIR0
DTDZBS0=IF_FS_DIR0

!similarly for C (ltdir)
DLDZCA=(R%NTMAX+2)/2
DLDZCS=(R%NTMAX+3)/2
DTDZCA0=IF_FS_DIR0
DTDZCS0=IF_FS_DIR0

! competition: NPRTRV ... larger == NUMP ... larger == NSMAX/NPRTRW
!              setting NPRTRV=20 ... leads to 7GB ZAA since NUMP==55

!allocate matrices for matrix multiplications
!ALLOCATE(IZBA(IF_FS_INV0*TDZAA*D%NUMP))
ALLOCATE(IZBS(IF_FS_INV0*TDZAS*D%NUMP))
print*,"New: allocating IZBS as a 1D array!"
! just use IZBS
!IZBA=>IZBS(:,1:TDZAA,:)

ALLOCATE(ZAA(R%NDGNH,TDZAA,D%NUMP))
ALLOCATE(ZAS(R%NDGNH,TDZAS,D%NUMP))

! Allocate matrices for rescaling to allow half-precision Legendre transforms
!ALLOCATE(ZAMAX(IF_FS_INV0,D%NUMP))
!ALLOCATE(ZSMAX(IF_FS_INV0,D%NUMP))

! transpose of C (for better memory access patterns)
!ALLOCATE(IZCAT(IF_FS_INV0,R%NDGNH,D%NUMP))
ALLOCATE(IZCST(IF_FS_INV0*R%NDGNH*D%NUMP))

!ALLOCATE(DZBAT(IF_FS_DIR0,R%NDGNH,D%NUMP))
ALLOCATE(DZBST(IF_FS_DIR0*R%NDGNH*D%NUMP))

! transpose of C (for better memory access patterns)
ALLOCATE(DZCAT(IF_FS_DIR0*TDZAA*D%NUMP))
ALLOCATE(DZCST(IF_FS_DIR0*TDZAS*D%NUMP))
DZCAT(:) = 0
DZCST(:) = 0
IZCST(:) = 0
!DZCAT=>DZCST(:,1:TDZAA,:)

write(nout,*)'sizes NUMP=',D%NUMP
write(nout,*)'ZAS:',size(ZAS)
write(nout,*)'IZBS :',size(IZBS )
write(nout,*)'IZCST:',size(IZCST)
write(nout,*)'DZBST:',size(DZBST)
write(nout,*)'DZCST:',size(DZCST)
write(nout,*)'DZCAT:',size(DZCAT)
!!!$ACC ENTER DATA CREATE(ZAA,ZAS,IZBA,IZBS,IZCAT,IZCST,DZBAT,DZBST,DZCAT,DZCST) 

write(nout,*) "Copy arrays and structures to GPU"
!$ACC ENTER DATA COPYIN(ZAA,ZAS,IZBS,IZCST,DZBST,DZCST,DZCAT) &
!$ACC& COPYIN(F,F%RN,F%RLAPIN,S,S%FA,S%ITHRESHOLD,S%LUSEFLT,D,D%NUMP,D%MYMS,R,R%NDGNH,R%NSMAX,G,G%NDGLU) &
!$ACC& copyin(D%NPNTGTB0,D%NPNTGTB1,D%NSTAGT0B,D%NSTAGT1B,D%NSTAGTF,G%NMEN,D%NPROCM,D%NPTRLS,G,G%NLOEN,D%MSTABF)

! Initialize A arrays

izbs = 0._JPRBT
!$acc update device(izbs)
dzbst = 0._JPRBT
!$acc update device(dzbst)

! zero arrays
!$ACC PARALLEL LOOP
DO JMLOC=1,D%NUMP
  !$ACC loop
  DO JK=1,TDZAA
    !$ACC loop
    DO J=1,LDZAA
      ZAA(J,JK,JMLOC)=0._JPRBT
    ENDDO
  ENDDO
ENDDO

!$ACC PARALLEL LOOP
DO JMLOC=1,D%NUMP
  !$ACC loop
  DO JK=1,TDZAS
    !$ACC LOOP
    DO J=1,LDZAS
      ZAS(J,JK,JMLOC)=0._JPRBT
    ENDDO
  ENDDO
ENDDO

! Do this on the host

zaa(:,:,:) = 0
DO JMLOC=1,D%NUMP
  KM = D%MYMS(JMLOC)   
  KDGLU = MIN(R%NDGNH,G%NDGLU(KM))
   
  ILA = (R%NSMAX-KM+2)/2
  DO JK=1,KDGLU
    DO J=1,ILA
      ZAA(JK,J,JMLOC)=S%FA(JMLOC)%RPNMA(JK,J)
    ENDDO
  ENDDO
ENDDO

ZAS(:,:,:) = 0
DO JMLOC=1,D%NUMP
  KM = D%MYMS(JMLOC)
  KDGLU = MIN(R%NDGNH,G%NDGLU(KM))

  ILS = (R%NSMAX-KM+3)/2
  DO JK=1,KDGLU
    DO J=1,ILS
      ZAS(JK,J,JMLOC)=S%FA(JMLOC)%RPNMS(JK,J)
    ENDDO
  ENDDO
ENDDO

! permanent copy of Legendre polynomials into device

!$ACC update device(ZAA)
!$ACC update device(ZAS)

IF_FOUBUF=MAX(IF_OUT_LT,IF_FS)
ALLOCATE(FOUBUF_IN(MAX(1,D%NLENGT0B*2*IF_FOUBUF)))
ALLOCATE(FOUBUF(MAX(1,D%NLENGT0B*2*IF_FOUBUF)))
! memory save

ALLOCATE(ZGTF(2*IF_FS,D%NLENGTF))
write(nout,*) 'Create ZGTF and Leg. arrays on GPU :',size(ZGTF)
!$ACC enter data create(ZGTF)

ALLOCATE(ZIA(IF_FS_INV0,R%NLEI1,D%NUMP))
ALLOCATE(ZEPSNM(d%nump,0:R%NTMAX+2))
ALLOCATE(ZSOA1(2*IF_OUT_LT,R%NLEI3,D%NUMP))
ALLOCATE(ZAOA1(2*IF_OUT_LT,R%NLEI3,D%NUMP))
ALLOCATE(ISTAN(D%NUMP,R%NDGNH))
ALLOCATE(ISTAS(D%NUMP,R%NDGNH))
!ALLOCATE(ZSIA(IF_FS_INV0,R%NDGNH,D%NUMP))
ALLOCATE(ZAIA(IF_FS_INV0,R%NDGNH,D%NUMP))
ALLOCATE(ZOA1(4*IF_FS_DIR0,R%NLED4,D%NUMP))
ALLOCATE(ZOA2(MAX(4*IF_UV,1),R%NLED4,D%NUMP))
write(nout,*)'ZIA  :',size(ZIA  )
write(nout,*)'ZSOA1:',size(ZSOA1)
write(nout,*)'ZAOA1:',size(ZAOA1)
write(nout,*)'ZAIA :',size(ZAIA )
write(nout,*)'ZOA1 :',size(ZOA1 )
write(nout,*)'ZOA2 :',size(ZOA2 )
!!!$ACC enter data create(ZIA,ZEPSNM,ZSOA1,ZAOA1,ISTAN,ISTAS,ZSIA,ZAIA,ZOA1,ZOA2)
!$ACC enter data create(ZIA,ZEPSNM,ZSOA1,ZAOA1,ZAIA,ZOA1,ZOA2)

zepsnm = 0._JPRBT
CALL PREPSNM
!$acc update device(zepsnm)
zgtf = 0._JPRBT
!$acc update device(zgtf)
zia = 0._JPRBT
!$acc update device(zia)
!zsia = 0._JPRBT
!!!$acc update device(zsia)
zaia = 0._JPRBT
!$acc update device(zaia)
zoa1 = 0._JPRBT
!$acc update device(zoa1)
zoa2 = 0._JPRBT
!$acc update device(zoa2)
zaoa1 = 0._JPRBT
!$acc update device(zaoa1)
zsoa1 = 0._JPRBT
!$acc update device(zsoa1)

! add arrays for GPNORM1
ALLOCATE(ZAVE(IF_FS,R%NDGL))
ALLOCATE(ZMINGL(IF_FS,R%NDGL))
ALLOCATE(ZMAXGL(IF_FS,R%NDGL))
ALLOCATE(ZMINGPN(IF_FS))
ALLOCATE(ZMAXGPN(IF_FS))
!$ACC enter data create(ZAVE,ZMINGL,ZMAXGL,ZMINGPN,ZMAXGPN)

zave = 0._JPRBT
!$acc update device(zave)
zmingl = 0._JPRBT
!$acc update device(zmingl)
zmaxgl = 0._JPRBT
!$acc update device(zmaxgl)
zmingpn = 0._JPRBT
!$acc update device(zmingpn)
zmaxgpn = 0._JPRBT
!$acc update device(zmaxgpn)

!set up flat copies of constant data
R_NSMAX=R%NSMAX
R_NTMAX=R%NTMAX
R_NDGNH=R%NDGNH
R_NDGL=R%NDGL
R_NNOEXTZL=R%NNOEXTZL


ALLOCATE(D_NSTAGT0B(SIZE(D%NSTAGT0B)))
ALLOCATE(D_NSTAGT1B(SIZE(D%NSTAGT1B)))
ALLOCATE(D_NPNTGTB0(0:SIZE(D%NPNTGTB0,1)-1,SIZE(D%NPNTGTB0,2)))
ALLOCATE(D_NPNTGTB1(SIZE(D%NPNTGTB1,1),SIZE(D%NPNTGTB1,2)))
ALLOCATE(D_MYMS(SIZE(D%MYMS)))
ALLOCATE(D_NPROCL(SIZE(D%NPROCL)))
ALLOCATE(D_NASM0(0:SIZE(D%NASM0)-1))
ALLOCATE(D_NSTAGTF(SIZE(D%NSTAGTF)))
ALLOCATE(D_MSTABF(SIZE(D%MSTABF)))
ALLOCATE(D_NPROCM(0:SIZE(D%NPROCM)-1))
ALLOCATE(D_NPTRLS(SIZE(D%NPTRLS)))

ALLOCATE(G_NDGLU(0:SIZE(G%NDGLU)-1))
ALLOCATE(G_NMEN(SIZE(G%NMEN)))
ALLOCATE(G_NLOEN(SIZE(G%NLOEN)))

ALLOCATE(F_RW(SIZE(F%RW)))


DO I=0,SIZE(G%NDGLU)-1
   G_NDGLU(I)=G%NDGLU(I)
end DO

G_NMEN_MAX=0
DO I=1,SIZE(G%NMEN)
   G_NMEN(I)=G%NMEN(I)
   if (G_NMEN(I) .gt. G_NMEN_MAX) G_NMEN_MAX=G_NMEN(I)
end DO

G_NLOEN_MAX=0
DO I=1,SIZE(G%NLOEN)
   G_NLOEN(I)=G%NLOEN(I)
   if (G_NLOEN(I) .gt. G_NLOEN_MAX) G_NLOEN_MAX=G_NLOEN(I)
end DO

DO I=1,SIZE(D%NSTAGT0B)
   D_NSTAGT0B(I)=D%NSTAGT0B(I)
END DO

DO I=1,SIZE(D%NSTAGT1B)
   D_NSTAGT1B(I)=D%NSTAGT1B(I)
END DO

DO I=1,SIZE(D%NPROCL)
   D_NPROCL(I)=D%NPROCL(I)
END DO

DO I=0,SIZE(D%NASM0)-1
   D_NASM0(I)=D%NASM0(I)
END DO

DO I=1,SIZE(D%NSTAGTF)
   D_NSTAGTF(I)=D%NSTAGTF(I)
END DO

DO I=1,SIZE(D%MSTABF)
   D_MSTABF(I)=D%MSTABF(I)
END DO

DO I=0,SIZE(D%NPROCM)-1
   D_NPROCM(I)=D%NPROCM(I)
END DO

DO I=1,SIZE(D%NPTRLS)
   D_NPTRLS(I)=D%NPTRLS(I)
END DO

DO I=1,SIZE(D%NPNTGTB0,2)
   DO J=0,SIZE(D%NPNTGTB0,1)-1
      D_NPNTGTB0(J,I)=D%NPNTGTB0(J,I)
   end DO
END DO

DO I=1,SIZE(D%NPNTGTB1,2)
   DO J=1,SIZE(D%NPNTGTB1,1)
      D_NPNTGTB1(J,I)=D%NPNTGTB1(J,I)
   end DO
END DO

D_NUMP=D%NUMP

KMLOC0 = -1
DO I=1,SIZE(D%MYMS)
   D_MYMS(I)=D%MYMS(I)
   IF(D_MYMS(I) == 0) KMLOC0 = I
end DO

! arrays for m=0 in ledir_mod:
IF(KMLOC0 >= 0) THEN
  ALLOCATE(ZAA0(R%NDGNH,TDZAA))
  ALLOCATE(ZAS0(R%NDGNH,TDZAS))
  ALLOCATE(DZBST0(IF_FS_DIR0*R%NDGNH))
  ALLOCATE(DZCAT0(IF_FS_DIR0*TDZAA))
  ALLOCATE(DZCST0(IF_FS_DIR0*TDZAS))
  DZCAT0(:) = 0
  DZCST0(:) = 0
  !$ACC ENTER DATA COPYIN(ZAA0,DZBST0,DZCAT0,ZAS0,DZCST0)
  ZAA0 = ZAA(:,:,KMLOC0)
  ZAS0 = ZAS(:,:,KMLOC0)
  !$ACC update device(ZAA0)
  !$ACC update device(ZAS0)
  dzbst0 = 0._JPRD
  !$acc update device(dzbst0)
  WRITE(NOUT,*) 'GPU arrays for m=0 successfully allocated'
ENDIF

DO I=1,SIZE(F%RW)
   F_RW(I)=F%RW(I)
END DO

!$ACC ENTER DATA COPYIN(R_NSMAX,R_NTMAX,R_NDGL,R_NNOEXTZL,R_NDGNH,D_NSTAGT0B,D_NSTAGT1B,&
!$ACC&                  D_NPNTGTB1,D_NPROCL,D_NUMP,D_MYMS,D_NASM0,D_NSTAGTF,D_MSTABF,&
!$ACC&                  D_NPNTGTB0,D_NPROCM,D_NPTRLS,G_NDGLU,G_NMEN,G_NMEN_MAX,G_NLOEN,&
!$ACC&                  G_NLOEN_MAX,F_RW)

WRITE(NOUT,*) 'GPU copy arrays (r, d, g) successfully allocated'
!$ACC wait

! free memory
!DO JMLOC=1,D%NUMP
!  DEALLOCATE(S%FA(JMLOC)%RPNMA)
!  DEALLOCATE(S%FA(JMLOC)%RPNMS)
!ENDDO

!endif INTERFACE

ENDIF

END SUBROUTINE SETUP_TRANS
