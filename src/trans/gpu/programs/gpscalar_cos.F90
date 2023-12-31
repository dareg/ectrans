! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

PROGRAM GPSCALAR_COS
! Nils Wedi 2010, ECMWF, test transform of a single scalar field to lat-lon

USE PARKIND_ECTRANS, ONLY : JPIM, JPRBT, JPRD, JPIB
USE YOMHOOK, ONLY : LHOOK
USE GRIB_API
USE MPL_MODULE
USE, INTRINSIC :: ISO_C_BINDING, ONLY:  C_PTR, C_INT, C_NULL_PTR,C_SIZE_T,C_F_POINTER
USE SHAREDMEM_MOD
USE BYTES_IO_MOD

IMPLICIT NONE

CHARACTER(LEN=13) :: OPTIONS
DATA OPTIONS/'g:G:l:t:sLOh;'/
CHARACTER(LEN=127) :: CLARG , OPTLET
CHARACTER(LEN=127) :: COUTGPF,  CINGPF, COUTSPF, CINSPEC,  CINTEMP, cltypeOfGrid
INTEGER(KIND=JPIM) ::  IFILE, IFILESP, IFILET, IGRIB_SP,IGRIB_H, IGRIB_OUT, IGRIB_T,ILEGFILE

INTEGER(KIND=JPIM) :: ISTACK,IMAX_THREADS
INTEGER(KIND=JPIM) :: getstackusage
INTEGER(KIND=JPIM) :: OPTVAL,OMP_GET_MAX_THREADS

INTEGER(KIND=JPIM) :: NERR,NSMAX,NDGLI, NDGLO, NSMAX_I, IDGL, NSIZE, ICOUNT
INTEGER(KIND=JPIM) :: NOUT,NSPEC2,NGPTOT,NGPTOTG,IMAXFLD,IFLD, NSPEC2G, ISTRUNC
INTEGER(KIND=JPIM) :: ITAG,IRET, IOUT, ILAT, ILON
INTEGER(KIND=JPIM) :: IOUTGPF, I
INTEGER(KIND=JPIM) ,ALLOCATABLE :: NLOEN(:),ITO(:),NPRCIDS(:)
INTEGER(KIND=JPIM) ,POINTER :: IIBUF(:)
INTEGER(KIND=JPIM) :: GETOPT,JJ
REAL(KIND=JPRBT),ALLOCATABLE :: ZSPEC(:,:),ZFPDAT(:),ZSPECG(:,:),ZSPECM(:,:)
REAL(KIND=JPRBT),ALLOCATABLE :: ZG(:,:,:), ZGG(:,:), ZGGOUT(:),ZGGM(:,:,:)
LOGICAL :: LLSCAL, LLATLON_OUT,LLATLON_IN, LLSPECIN,LLLL

INTEGER(KIND=JPIM) :: IPOI, JL, JGL, J, IR
REAL(KIND=JPRBT)    :: ZDEG2RAD, ZLON, ZLONI
REAL(KIND=JPRBT)    :: RPI, RA
TYPE(C_PTR)        :: MEMBUF
INTEGER(C_SIZE_T)     :: IBYTES,IOFF
INTEGER(KIND=JPIM) :: IPROC, NPROC, MYPROC, IMYSETV,IRESOL
INTEGER(KIND=JPIM) :: ILONS, IPLPRESENT, ITEST
REAL(KIND=JPRD)    :: ZTIME0,ZTIME1,ZTIME2


!     ------------------------------------------------------------------

#include "setup_trans0.h"
#include "setup_trans.h"
#include "inv_trans.h"
#include "dir_trans.h"
#include "dist_grid.h"
#include "dist_spec.h"
#include "gath_grid.h"
#include "gath_spec.h"
#include "trans_inq.h"
#include "trans_release.h"
#include "user_clock.h"

! Initializations
LHOOK=.TRUE.
NERR = 0
NOUT = 6
IMAXFLD = 1

! Set defaults for options

CINSPEC = 'insp.grib'
CINGPF  = 'in.grib'
CINTEMP = 'template.grib'
COUTGPF = 'out.grib'
COUTSPF = 'outspec.grib'
IPROC   = 1
LLSCAL   = .false.
NDGLO    = 0
NSMAX = 0
LLATLON_OUT=.false.
LLATLON_IN=.false.
LLSPECIN=.false.

RPI = 2.0_JPRBT*ASIN(1.0_JPRBT)
RA  = 6371229._JPRBT

! Crack options

DO
  OPTVAL = GETOPT(OPTIONS,CLARG)
  IF(OPTVAL <= 0) EXIT
  OPTLET=CHAR(OPTVAL)
  IF(OPTVAL <= 0) EXIT
  IF(OPTLET /= 'g'.AND.&
   &OPTLET  /= 'l'.AND.&
   &OPTLET  /= 't'.AND.&
   &OPTLET  /= 's'.AND.&
   &OPTLET  /= 'h'.AND.&
   &OPTLET  /= 'L'.AND.&
   &OPTLET  /= 'O'.AND.&
   &OPTLET  /= 'G')  THEN
    CALL USE_GPT
    CALL ABOR1('GPSCAL_DER:ERROR IN ARGUMENT')
  ENDIF
  IF(OPTLET == 'h') THEN
    CALL USE_GPT
    STOP
  ELSEIF(OPTLET == 'l')  THEN
    READ(CLARG,*) NDGLO
  ELSEIF(OPTLET == 't')  THEN
    READ(CLARG,*) NSMAX
  ELSEIF(OPTLET == 'g')  THEN
    CINGPF = CLARG
  ELSEIF(OPTLET == 'G') THEN
    COUTGPF = CLARG
  ELSEIF(OPTLET == 's') THEN
    LLSPECIN=.true.
  ELSEIF(OPTLET == 'L') THEN
    LLATLON_IN=.true.
  ELSEIF(OPTLET == 'O') THEN
    LLATLON_OUT=.true.
  ENDIF
ENDDO

IF( NDGLO == 0 ) THEN
  WRITE(NERR,*) 'NDGLO =',NDGLO
  CALL ABOR1('GPT: NUMBER OF GAUSSIAN LATITUDES NOT SET (use -l option)')
ENDIF
IF( NSMAX == 0 ) THEN
  WRITE(NERR,*) 'NSMAX =',NSMAX
  CALL ABOR1('GPT: TRUNCATION NOT SET (use -t option)')
ENDIF

! Message passing setup
! Participating processors limited by -P option
IF(IPROC > 1 ) THEN
  CALL MPL_INIT
  NPROC   = MPL_NPROC()
  MYPROC = MPL_MYRANK()
ELSE
  NPROC  = 1
  MYPROC = 1
  LLSCAL=.true.
ENDIF

ALLOCATE(NPRCIDS(NPROC))
DO JJ=1,NPROC
  NPRCIDS(JJ) = JJ
ENDDO

!-------------------------
ITAG = 191919
IMAX_THREADS = -99
!$ IMAX_THREADS = OMP_GET_MAX_THREADS()
print *,' IMAX_THREADS = ',IMAX_THREADS
!------------------------------------------------------
! INPUT GRIDPOINT FIELD JUST TO SETUP OUTPUT GG or LAT-LON
!------------------------------------------------------
! open file
CALL GRIB_OPEN_FILE(IFILE,CINGPF,'R')
! get handle
CALL GRIB_NEW_FROM_FILE(IFILE,IGRIB_H, IRET)
! get key and check
CALL GRIB_GET(IGRIB_H,'typeOfGrid', cltypeOfGrid)
!!IF ( TRIM(cltypeOfGrid) /='reduced_gg' ) THEN
!  CALL ABOR1('GPT: Need gridpoint field in output resolution as input!')
!ENDIF
CALL GRIB_GET(IGRIB_H,'Nj',NDGLI)
CALL GRIB_GET_SIZE(IGRIB_H,'values',NSIZE)
!CALL GRIB_GET(KGRIB_HANDLE,'latitudeOfFirstGridPoint',INP)
!CALL GRIB_GET(KGRIB_HANDLE,'latitudeOfLastGridPoint',ISP)
!ZDLAT = ZDTOR*REAL(ISP-INP,JPRBT)/REAL(ILATS-1,JPRBT)
ALLOCATE(NLOEN(NDGLI))
CALL GRIB_GET(IGRIB_H,'PLPresent',IPLPresent)
IF (IPLPresent == 0) THEN
  CALL GRIB_GET(IGRIB_H,'numberOfPointsAlongAParallel',ILONS)
  NLOEN(:) = ILONS
ELSE
  CALL GRIB_GET(IGRIB_H,'pl',NLOEN)
ENDIF

!------------------------------------------------------
! TEMPLATE GRIDPOINT FIELD JUST TO SETUP OUTPUT (could link to in.grib)
!------------------------------------------------------
CALL GRIB_OPEN_FILE(IFILET,CINTEMP,'R')
CALL GRIB_NEW_FROM_FILE(IFILET,IGRIB_T, IRET)
CALL GRIB_CLOSE_FILE(IFILET)  

!------------------------------------------------------
! INPUT SPECTRAL FIELD IS OPTIONAL AND WILL THEN JUST PERFORM INV_TRANS
!------------------------------------------------------

IF( LLSPECIN ) THEN
  ! open file
  CALL GRIB_OPEN_FILE(IFILESP,CINSPEC,'R')
  ! get handle
  CALL GRIB_NEW_FROM_FILE(IFILESP,IGRIB_SP,IRET)
  ! get key and check
  CALL GRIB_GET(IGRIB_SP,'typeOfGrid', cltypeOfGrid)
  IF ( TRIM(cltypeOfGrid) /= 'sh' ) THEN
    CALL ABOR1('GPT: Need spectral field on input!')
  ENDIF
  CALL GRIB_GET(IGRIB_SP,'J',NSMAX_I)
  IF( NSMAX_I /= NSMAX ) THEN
    CALL ABOR1('GPT: NSMAX in file different from specified!')
  ENDIF
ENDIF

!------------------------------------------------------

! Prepare for transforms
CALL USER_CLOCK(PELAPSED_TIME=ZTIME0)
write(0,*) ' BEFORE SETUP_TRANS0' ,ZTIME0
CALL SETUP_TRANS0(KOUT=NOUT,KERR=0,KPRINTLEV=0,KMAX_RESOL=1,&
 &                KPRGPNS=NPROC,KPRGPEW=1,KPRTRW=NPROC,LDMPOFF=LLSCAL)

! setup regular Gaussian grid with lat/lon equivalent dual
IF( NDGLI > NDGLO+1 ) THEN
  CALL ABOR1('GPT: NDGLI > NDGLO+1 not possible !')
ENDIF
LLLL = LLATLON_IN .OR. LLATLON_OUT
IBYTES=3291640360_C_SIZE_T
CALL SHAREDMEM_MALLOC_BYTES(MEMBUF,IBYTES)
CALL BYTES_IO_OPEN(ILEGFILE,'legpol_T1279','R')
CALL C_F_POINTER(MEMBUF,IIBUF,(/IBYTES/4/))
CALL BYTES_IO_READ(ILEGFILE,IIBUF,IBYTES,IRET)
CALL USER_CLOCK(PELAPSED_TIME=ZTIME0)
write(0,*) ' AFTER READING LEGPOL' ,ZTIME0
!IOFF = 1
!DO 
!  CALL BYTES_IO_READ(ILEGFILE,IIBUF(IOFF:),1024*1024,IRET)
!  WRITE(0,*) 'IOFF= ',IOFF
!  IOFF=IOFF+1024*1024/4
!  IF(IRET == GRIB_END_OF_FILE) EXIT
!ENDDO
CALL BYTES_IO_CLOSE(ILEGFILE,IRET)
CALL USER_CLOCK(PELAPSED_TIME=ZTIME0)
write(0,*) ' BEFORE SETUP_TRANS' ,ZTIME0
IF( NDGLI == NDGLO+1 ) THEN
  IDGL=NDGLO+2
  CALL SETUP_TRANS(KRESOL=IRESOL,KSMAX=NSMAX,KDGL=NDGLO,LDSPLIT=.FALSE.,KLOEN=NLOEN,LDUSEFLT=.false.,&
!   & LDLL=LLLL,CDIO_LEGPOL='membuf',KLEGPOLPTR=MEMBUF,KLEGPOLPTR_LEN=IBYTES)
   & LDLL=LLLL,CDIO_LEGPOL='writef',CDLEGPOLFNAME='legpol_T159_bf')
!   & LDLL=LLLL,CDIO_LEGPOL='readf',CDLEGPOLFNAME='legpol_T1279')
!   & LDLL=LLLL)
ELSE
  IDGL=NDGLO
  CALL SETUP_TRANS(KRESOL=IRESOL,KSMAX=NSMAX,KDGL=NDGLO,LDSPLIT=.FALSE.,KLOEN=NLOEN,LDUSEFLT=.false.,&
!   & LDLL=LLLL,LDSHIFTLL=LLLL,CDIO_LEGPOL='membuf',KLEGPOLPTR=MEMBUF,KLEGPOLPTR_LEN=IBYTES)
!   & LDLL=LLLL,LDSHIFTLL=LLLL,CDIO_LEGPOL='readf',CDLEGPOLFNAME='legpol_T159_bf')
   & LDLL=LLLL,LDSHIFTLL=LLLL,CDIO_LEGPOL='writef',CDLEGPOLFNAME='legpol_T159_bf')
!   & LDLL=LLLL,LDSHIFTLL=LLLL)
ENDIF
CALL USER_CLOCK(PELAPSED_TIME=ZTIME0)
write(0,*) ' AFTER SETUP_TRANS' ,ZTIME0

CALL TRANS_INQ(KRESOL=IRESOL,KSPEC2=NSPEC2,KSPEC2G=NSPEC2G,KGPTOT=NGPTOT,KGPTOTG=NGPTOTG,KMYSETV=IMYSETV)

write(NERR,*) 'LATITUDES INPUT, LATITUDES OUPUT, RESOL', NDGLI, NDGLO, NSMAX
write(NERR,*) 'DIMS ', NSPEC2,NSPEC2G,NGPTOT,NGPTOTG

ALLOCATE(ZSPEC(IMAXFLD,NSPEC2))
ALLOCATE(ZG(NGPTOT,1,1))

IF( MYPROC == 1 ) THEN
  ALLOCATE(ZSPECG(1,NSPEC2G))
  ALLOCATE(ZGG(NGPTOTG,IMAXFLD))
ELSE
  ALLOCATE(ZSPECG(1,NSPEC2))
  ALLOCATE(ZGG(NGPTOT,IMAXFLD))
ENDIF

ALLOCATE(ITO(IMAXFLD))

! Gridpoint to spectral transform

ITO(:) = 1
IFLD = 0

ITEST = 0
IFLD=1

IF( .NOT.LLSPECIN ) THEN
  IF( MYPROC == 1 ) THEN
    IF( ITEST == 1 ) THEN
      ALLOCATE(ZFPDAT(NGPTOTG))
      ! test function COS(lambda)
      IPOI=0
      ZDEG2RAD=RPI/180._JPRBT
      DO JGL=1,NDGLI
        ZLONI=360._JPRBT/FLOAT(NLOEN(JGL))*ZDEG2RAD
        DO JL=1,NLOEN(JGL)
          IPOI=IPOI+1
          ZLON=FLOAT(JL-1)*ZLONI
          ZFPDAT(IPOI)=COS(ZLON)
        ENDDO
      ENDDO
      ZGG(:,1) = ZFPDAT(:)
      DEALLOCATE(ZFPDAT) 
    ELSE

      !needs mapping to field suitable for transforms
      IF( NDGLI == NDGLO+1 ) THEN
        ! odd number of latitudes, lat/lon field
        ALLOCATE(ZGGOUT(NSIZE))
        ! read input data
        CALL GRIB_GET(IGRIB_H,'values',ZGGOUT)
        ICOUNT=0
        IR=0
        DO ILAT=1,IDGL
          DO ILON=1,2*NDGLO
            ! need equator 2x in input data, duplicate one row
            IF( ILAT==IDGL/2+1 .AND. IR==0 ) THEN
              ICOUNT=ICOUNT-2*NDGLO
              IR=1
            ENDIF
            ICOUNT=ICOUNT+1
            ZGG(ILON+(ILAT-1)*2*NDGLO,1) =  ZGGOUT(ICOUNT)
          ENDDO
        ENDDO
        DEALLOCATE(ZGGOUT)
      ELSE
        !  standard dimensions
        CALL GRIB_GET(IGRIB_H,'values',ZGG(1:NGPTOTG,1))
      ENDIF
    ENDIF
  ENDIF

  ! Distribute gridpoint fields to processors
  CALL DIST_GRID(PGPG=ZGG,KFDISTG=IFLD,KFROM=ITO,PGP=ZG)
  
  ! Spectral transform
  CALL DIR_TRANS(KRESOL=IRESOL,PSPSCALAR=ZSPEC(1:IFLD,:),PGP=ZG,LDLATLON=LLATLON_IN)
  
  ! Gather spectral fields to processor 1
  CALL GATH_SPEC(KRESOL=IRESOL,PSPECG=ZSPECG,KFGATHG=1,KTO=ITO,PSPEC=ZSPEC)

  IF(MYPROC == 1) THEN
    !  -----------------------------
    
    !*    WRITE THE SPECTRAL FIELD
    
    ALLOCATE(ZFPDAT(NSPEC2G))
    CALL GRIB_NEW_FROM_SAMPLES(IGRIB_OUT,'sh_sfc_grib2')
    CALL GRIB_OPEN_FILE(IOUT,COUTSPF,'w')
    CALL GRIB_SET(IGRIB_OUT,'gridType','sh')
    CALL GRIB_SET(IGRIB_OUT,'numberOfBitsContainingEachPackedValue',16)
    CALL GRIB_SET(IGRIB_OUT,'pentagonalResolutionParameterJ',NSMAX)
    CALL GRIB_SET(IGRIB_OUT,'pentagonalResolutionParameterK',NSMAX)
    CALL GRIB_SET(IGRIB_OUT,'pentagonalResolutionParameterM',NSMAX)
    CALL GRIB_SET(IGRIB_OUT,'laplacianOperator',0.5_JPRBT)
    IF(NSMAX>= 213) THEN
      ISTRUNC = 20
    ELSE
      ISTRUNC = MIN(10,NSMAX)
    ENDIF
    CALL GRIB_SET(IGRIB_OUT,'subSetJ',ISTRUNC)
    CALL GRIB_SET(IGRIB_OUT,'subSetK',ISTRUNC)
    CALL GRIB_SET(IGRIB_OUT,'subSetM',ISTRUNC)
    CALL GRIB_SET(IGRIB_OUT,'paramId',80)
    CALL GRIB_SET(IGRIB_OUT,'level',0)
    CALL GRIB_SET(IGRIB_OUT,'values',ZSPECG(1,:))
    CALL GRIB_WRITE(IGRIB_OUT,IOUT)
    CALL GRIB_RELEASE(IGRIB_OUT)
    CALL GRIB_CLOSE_FILE(IOUT)
    DEALLOCATE(ZFPDAT)
    
  ENDIF
  
  ! Syncronize processors
  IF(NPROC > 1) THEN
    CALL MPL_BARRIER(CDSTRING='GPSCALAR_COS:')
  ENDIF

ELSE
  IF(MYPROC == 1) THEN
    ! read input data
    CALL GRIB_GET(IGRIB_SP,'values',ZSPECG(1,:))
  ENDIF
  CALL DIST_SPEC(PSPECG=ZSPECG,KFDISTG=IFLD,KFROM=ITO,PSPEC=ZSPEC)
ENDIF

!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     LDSCDERS    - indicating if derivatives of scalar variables are req.
!     PGP(:,:,:) - gridpoint fields (output)
!                  PGP need to  dimensioned (NPROMA,IF_GP,NGPBLKS) where
!                  NPROMA is the blocking factor, IF_GP the total number
!                  of output fields and NGPBLKS the number of NPROMA blocks.
!                  The ordering of the output fields is as follows (all 
!                  parts are optional depending on the input switches):
!       scalar fields : IF_SCALARS_G fields (if pspscalar present)
!       N-S derivative of scalar fields : IF_SCALARS_G fields (if pspscalar
!                                         present and LDSCDERS)
!       E-W derivative of scalar fields : IF_SCALARS_G fields (if pspscalar
!                                         present and LDSCDERS)
!       IF_SCALARS_G is the GLOBAL number of scalar fields as giben by the 
!       length of KVESETSC (or by number of fields in PSPSCALAR if no spectral
!       'b-set' split

! inverse transform

!do jgl=1,100
CALL INV_TRANS(KRESOL=IRESOL,PSPSCALAR=ZSPEC(1:IFLD,:),PGP=ZG,LDLATLON=LLATLON_OUT)
!enddo
CALL USER_CLOCK(PELAPSED_TIME=ZTIME0)
write(0,*) ' AFTER INV_TRANS' ,ZTIME0
!ALLOCATE(ZG(NGPTOT,1,1))
!allocate(zggm(NGPTOT,100,1))
!ALLOCATE(ZSPECM(100,NSPEC2))
!ZSPECM=0
!CALL INV_TRANS(PSPSCALAR=ZSPECM,PGP=ZGGM,LDLATLON=LLATLON_OUT)
! Gather gridpoint fields to processor 1
CALL GATH_GRID(PGPG=ZGG,KFGATHG=IFLD,KTO=ITO,PGP=ZG)
CALL TRANS_RELEASE(KRESOL=IRESOL)
IF(MYPROC == 1) THEN

  !*    WRITE THE GRIDPOINT FIELD

  ! open file
  CALL GRIB_OPEN_FILE(IOUTGPF,COUTGPF,'W')
  print *,'IFLD ',IFLD
  CALL GRIB_CLONE(IGRIB_T,IGRIB_OUT)
  CALL GRIB_SET(IGRIB_OUT,'paramId',80)

  IF( NDGLI == NDGLO+1 ) THEN
    ! odd number of latitudes
    ALLOCATE(ZGGOUT(NSIZE))
    ICOUNT=0
    DO ILAT=1,IDGL
      DO ILON=1,2*NDGLO 
        ! equator is 2x in output data, remove one or both
        IF( ILAT/=IDGL/2+1 .AND. .NOT.(ILAT==IDGL/2.AND..NOT.LLATLON_OUT)) THEN
          ICOUNT=ICOUNT+1
          ZGGOUT(ICOUNT) = ZGG(ILON+(ILAT-1)*2*NDGLO,1)
        ENDIF
      ENDDO
    ENDDO
    write(NERR,*) 'OUTPUT DIMENSION: ', ICOUNT
    CALL GRIB_SET(IGRIB_OUT,'values',ZGGOUT(1: ICOUNT))
    DEALLOCATE(ZGGOUT) 
 ELSE
    !  standard dimensions
    CALL GRIB_SET(IGRIB_OUT,'values',ZGG(1:NGPTOTG,1))
  ENDIF
 

  CALL GRIB_WRITE(IGRIB_OUT,IOUTGPF)
  CALL GRIB_RELEASE(IGRIB_OUT)
  ! close output file
  CALL GRIB_CLOSE_FILE(IOUTGPF)  
ENDIF

! Syncronize processors
IF(NPROC > 1) THEN
  CALL MPL_BARRIER(CDSTRING='GPSCALAR_COS:')
ENDIF
IFLD = 0
  
!           gather stack usage statistics
ISTACK = GETSTACKUSAGE()
IF(MYPROC == 1) THEN
  PRINT 9000, istack
9000 FORMAT("Stack Utilisation Information",/,&
      &"=============================",//,&
      &"Node           Size(Bytes)",/,&
      &"====           ===========",//,&
      &"   1",11x,I10)
  
  DO I=2,NPROC
    CALL MPL_RECV(ISTACK,KSOURCE=NPRCIDS(I),KTAG=I, &
     & CDSTRING='GPSCALAR_COS:')
    PRINT '(I4,11X,I10)', I,ISTACK
  ENDDO
ELSE
  CALL MPL_SEND(ISTACK,KDEST=NPRCIDS(1),KTAG=MYPROC, &
   &   CDSTRING='GPSCALAR_COS:')
ENDIF

!Close down message passing
!--------------------------
CALL MPL_BARRIER(CDSTRING='GPSCALAR_COS:')
CALL MPL_END()
!--------------------------
STOP

CONTAINS  

!     ------------------------------------------------------------------

SUBROUTINE USE_GPT

!   WRITE MESSAGE TO INFORM ABOUT CORRECT USE OF PROGRAM
WRITE(NERR,*) ' CORRECT OPTIONS AND ARGUMENTS FOR GPSCALAR_DER ARE :'
WRITE(NERR,*) '   -g    input filename (grid-point) (def: fort.11)'
WRITE(NERR,*) '   -G    output grid-point filename (def: fort.10)'
WRITE(NERR,*) '   -l    ndgl'
WRITE(NERR,*) '   -h    displays options and arguments'
END SUBROUTINE USE_GPT

!     ------------------------------------------------------------------

END PROGRAM GPSCALAR_COS

