! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE SETUP_TRANS0(KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR,&
&                       KPRGPNS,KPRGPEW,KPRTRW,KCOMBFLEN,&
&                       LDMPOFF,LDSYNC_TRANS,KTRANS_SYNC_LEVEL,&
&                       LDEQ_REGIONS,K_REGIONS_NS,K_REGIONS_EW,K_REGIONS,&
&                       PRAD,LDALLOPERM,KOPT_MEMORY_TR)

!**** *SETUP_TRANS0* - General setup routine for transform package

!     Purpose.
!     --------
!     Resolution independent part of setup of transform package
!     Has to be called BEFORE SETUP_TRANS

!**   Interface.
!     ----------
!     CALL SETUP_TRANS0(...)

!     Explicit arguments : All arguments are optional, [..] default value
!     -------------------
!     KOUT - Unit number for listing output [6]
!     KERR - Unit number for error messages [0]
!     KPRINTLEV - level of output to KOUT, 0->no output,1->normal,2->debug [0]
!     KMAX_RESOL - maximum number of different resolutions for this run [1]
!     KPRGPNS - splitting level in N-S direction in grid-point space [1]
!     KPRGPEW - splitting level in E-W direction in grid-point space [1]
!     KPRTRW  - splitting level in wave direction in spectral space [1]
!     KCOMBFLEN - Size of communication buffer [1800000 (*8bytes) ]
!     LDMPOFF - switch off message passing [false]
!     LDSYNC_TRANS - switch to activate barriers in trmtol trltom [false]
!     KTRANS_SYNC_LEVEL - use of synchronization/blocking [0]
!     LDEQ_REGIONS - true if new eq_regions partitioning [false]
!     K_REGIONS    - Number of regions (1D or 2D partitioning)
!     K_REGIONS_NS - Maximum number of NS partitions
!     K_REGIONS_EW - Maximum number of EW partitions
!     PRAD         - Radius of the planet
!     LDALLOPERM  - Allocate certain arrays permanently
!     The total number of (MPI)-processors has to be equal to KPRGPNS*KPRGPEW

!     Method.
!     -------

!     Externals.  SUMP_TRANS0 - initial setup routine
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        R. El Khatib 03-01-24 LDMPOFF
!        G. Mozdzynski 2006-09-13 LDEQ_REGIONS
!        N. Wedi  2009-11-30 add radius

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NERR, NOUT, LMPOFF, LSYNC_TRANS, NTRANS_SYNC_LEVEL, MSETUP0, &
     &                      NMAX_RESOL, NPRINTLEV, NPROMATR, LALLOPERM,NSTACK_MEMORY_TR
USE TPM_DISTR       ,ONLY : LEQ_REGIONS, NCOMBFLEN, NPRGPEW,NPRGPNS, NPRTRW, NPRTRV, MYSETV
USE TPM_CONSTANTS   ,ONLY : RA
USE MPL_MODULE

USE SUMP_TRANS0_MOD ,ONLY : SUMP_TRANS0
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
USE EQ_REGIONS_MOD  ,ONLY : N_REGIONS, N_REGIONS_EW, N_REGIONS_NS
#ifdef _OPENACC
use openacc
#endif
use ec_env_mod, only : ec_getenv

! daand: added this to enable synchronization status
use cudafor

!endif INTERFACE

IMPLICIT NONE

INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN)  :: KOUT,KERR,KPRINTLEV,KMAX_RESOL,KPROMATR
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN)  :: KPRGPNS,KPRGPEW,KPRTRW,KCOMBFLEN
LOGICAL            ,OPTIONAL,INTENT(IN)  :: LDMPOFF
LOGICAL            ,OPTIONAL,INTENT(IN)  :: LDSYNC_TRANS
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN)  :: KTRANS_SYNC_LEVEL
LOGICAL            ,OPTIONAL,INTENT(IN)  :: LDEQ_REGIONS
LOGICAL            ,OPTIONAL,INTENT(IN)  :: LDALLOPERM
REAL(KIND=JPRB)    ,OPTIONAL,INTENT(IN)  :: PRAD
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(IN)  :: KOPT_MEMORY_TR
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT) :: K_REGIONS(:)
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT) :: K_REGIONS_NS
INTEGER(KIND=JPIM) ,OPTIONAL,INTENT(OUT) :: K_REGIONS_EW

INTEGER(KIND=JPIM) :: MYPROC
INTEGER :: IDEVICE_NUM,IDEVICE_TYPE, IPROC_PERNODE
integer :: idevtype, numdevs, mygpu, IERROR, istat
CHARACTER(LEN=2)  :: CL_NPROC_PERNODE

!ifndef INTERFACE

LOGICAL :: LLP1,LLP2

!     ------------------------------------------------------------------

MYPROC = MPL_MYRANK()


#ifdef gnarls

#ifdef _OPENACC
idevtype=acc_get_device_type()
numdevs = acc_get_num_devices(idevtype)
mygpu = mod(MYPROC-1,numdevs)
CALL acc_set_device_num(mygpu, idevtype)
mygpu = acc_get_device_num(idevtype)
write(*,*) 'MYPROC:',MYPROC, 'GPU:', mygpu, 'of ', numdevs
#endif

CL_NPROC_PERNODE=' '
CALL EC_GETENV('NPROC_PERNODE',CL_NPROC_PERNODE)
IF( CL_NPROC_PERNODE /= ' ')THEN
  READ(CL_NPROC_PERNODE,*) IPROC_PERNODE
  IDEVICE_NUM=MOD(MYPROC-1,IPROC_PERNODE)
  WRITE(0,'("TRANSFORM TEST: MYPROC=",I8," CL_NPROC_PERNODE=",A," IPROC_PERNODE=",I2,&
   & " IDEVICE_NUM=",I2)') MYPROC,CL_NPROC_PERNODE,IPROC_PERNODE,IDEVICE_NUM
  IDEVICE_TYPE=0
  !!CALL ACC_SET_DEVICE_NUM(IDEVICE_NUM,ACC_DEVICE_NVIDIA)
  CALL ACC_SET_DEVICE_NUM(IDEVICE_NUM,idevtype)
  !!CALL ACC_INIT(ACC_DEVICE_NVIDIA)
  CALL ACC_INIT(idevtype)
  !$OMP PARALLEL
  !!CALL ACC_SET_DEVICE_NUM(IDEVICE_NUM,ACC_DEVICE_NVIDIA)
  CALL ACC_SET_DEVICE_NUM(IDEVICE_NUM,idevtype)
  !!CALL ACC_INIT(ACC_DEVICE_NVIDIA)
  CALL ACC_INIT(idevtype)
!$OMP END PARALLEL
ENDIF
#endif



istat = cudaDeviceSynchronize()
write (*,*) "cudaDeviceSynchronize returned code ",istat


IF(MSETUP0 /= 0) THEN
!gr  CALL ABORT_TRANS('SETUP_TRANS0: SETUP_TRANS0 MAY ONLY BE CALLED ONCE')
ENDIF

! Default values

NOUT = 6
NERR = 0
NPRINTLEV = 0
NMAX_RESOL = 1
NPRGPNS = 1
NPRGPEW = 1
NPRTRW = 1
N_REGIONS_NS=1
N_REGIONS_EW=1
NPROMATR = 0
NCOMBFLEN = 1800000
LMPOFF = .FALSE.
LSYNC_TRANS=.FALSE.
NTRANS_SYNC_LEVEL=0
LEQ_REGIONS=.FALSE.
RA=6371229._JPRB
LALLOPERM=.FALSE.

! Optional arguments

IF(PRESENT(KOUT)) THEN
  NOUT = KOUT
ENDIF
IF(PRESENT(KERR)) THEN
  NERR = KERR
ENDIF
IF(PRESENT(KPRINTLEV)) THEN
  NPRINTLEV = KPRINTLEV
ENDIF

LLP1 = NPRINTLEV>0
LLP2 = NPRINTLEV>1
IF(LLP1) WRITE(NOUT,*) '=== ENTER ROUTINE SETUP_TRANS0 ==='

IF(PRESENT(KMAX_RESOL))THEN
  NMAX_RESOL = KMAX_RESOL
ENDIF
IF(PRESENT(KPROMATR))THEN
  IF(MOD(KPROMATR,2) /= 0) THEN
    CALL ABORT_TRANS('SETUP_TRANS0: KPROMATR HAS TO BE MULTIPLE OF 2')
  ENDIF
  NPROMATR = KPROMATR
ENDIF
IF(PRESENT(KPRGPNS)) THEN
  NPRGPNS = KPRGPNS
ENDIF
IF(PRESENT(KPRGPEW)) THEN
  NPRGPEW = KPRGPEW
ENDIF
IF(PRESENT(KPRTRW)) THEN
  NPRTRW = KPRTRW
ENDIF
IF(PRESENT(KCOMBFLEN)) THEN
  NCOMBFLEN = KCOMBFLEN
ENDIF
IF(PRESENT(LDMPOFF)) THEN
  LMPOFF = LDMPOFF
ENDIF
IF(PRESENT(LDSYNC_TRANS)) THEN
  LSYNC_TRANS = LDSYNC_TRANS
ENDIF
IF(PRESENT(KTRANS_SYNC_LEVEL)) THEN
  NTRANS_SYNC_LEVEL = KTRANS_SYNC_LEVEL
ENDIF
IF(PRESENT(LDEQ_REGIONS)) THEN
  LEQ_REGIONS = LDEQ_REGIONS
ENDIF
IF(PRESENT(KOPT_MEMORY_TR))THEN
  NSTACK_MEMORY_TR = KOPT_MEMORY_TR
ENDIF

! Initial setup
CALL SUMP_TRANS0

IF(PRESENT(K_REGIONS_NS)) THEN
  K_REGIONS_NS = N_REGIONS_NS
ENDIF

IF(PRESENT(K_REGIONS_EW)) THEN
  K_REGIONS_EW = N_REGIONS_EW
ENDIF

IF(PRESENT(K_REGIONS)) THEN
  IF(UBOUND(K_REGIONS,1) < N_REGIONS_NS) THEN
    write(0,*) "k_regions:",size(k_regions),":",k_regions
    CALL ABORT_TRANS('SETUP_TRANS0: K_REGIONS TOO SMALL')
  ELSE
    K_REGIONS(1:N_REGIONS_NS)=N_REGIONS(1:N_REGIONS_NS)
  ENDIF
ENDIF

IF(PRESENT(PRAD)) THEN
  RA=PRAD
ENDIF

IF(PRESENT(LDALLOPERM)) THEN
  LALLOPERM=LDALLOPERM
ENDIF

! Setup level 0 complete
MSETUP0 = 1

!     ------------------------------------------------------------------

!endif INTERFACE

END SUBROUTINE SETUP_TRANS0


