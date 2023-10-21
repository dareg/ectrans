! (C) Copyright 2000- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE FTDIR_MOD
CONTAINS
SUBROUTINE FTDIR(kf_fs)

!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFT992 - FFT routine
!     ----------
!
!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-04-24 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW

USE PARKIND_ECTRANS,ONLY : JPIM, JPIB, JPRBT
USE TPM_GEN,ONLY : NOUT
USE TPM_GEOMETRY,ONLY : G,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX
USE TPM_DIM,ONLY : R,R_NNOEXTZL
USE TPM_DISTR,ONLY : D,MYSETW,MYPROC,NPROC,D_NPTRLS,D_NSTAGTF,D_MSTABF,D_NSTAGT1B,&
  D_NPNTGTB0,D_NPROCM,D_NPROCL
USE TPM_TRANS,ONLY : ZGTF,FOUBUF_IN
USE TPM_FFTC,ONLY : CREATE_PLAN_FFT
USE CUDA_DEVICE_MOD

IMPLICIT NONE

integer(kind=jpim),intent(in) :: kf_fs

INTEGER(KIND=JPIM) :: KGL,IGLG,IST,ILEN,JM,JJ,JF,IOFF,inf,IPROC,ISTA
INTEGER(KIND=JPIM) :: IPLAN_R2C
INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,nlon
integer :: istat
real(kind=jprbt), allocatable :: zgtf2(:,:)

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

inf = size(zgtf,1)

!$ACC DATA COPYIN(IBEG,IEND,IINC,inf) PRESENT(ZGTF,FOUBUF_IN,D_NSTAGTF,D_NPTRLS,G_NMEN,&
!$ACC G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX,R_NNOEXTZL,D_NPROCM,D_NSTAGT1B,D_MSTABF,D_NPNTGTB0)

allocate(zgtf2(inf,size(zgtf,2)))

!$ACC DATA CREATE(ZGTF2)

DO KGL=IBEG,IEND,IINC
  IGLG=D_NPTRLS(MYSETW)+KGL-1
  IOFF=D_NSTAGTF(KGL)

  CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,G%NLOEN(IGLG),inf)

  !$ACC host_data use_device(ZGTF,ZGTF2)
  CALL EXECUTE_PLAN_FFTC(IPLAN_R2C,-1,ZGTF(1,IOFF+1),ZGTF2(1,IOFF+1))
  !$ACC end host_data
END DO

istat = cuda_Synchronize()

!$acc kernels DEFAULT(NONE)
zgtf(:,:) = zgtf2(:,:)
!$acc end kernels

!$acc end data

!$ACC parallel loop collapse(3) private(IGLG,nlon,IOFF,IST) DEFAULT(NONE)
DO KGL=IBEG,IEND,IINC
   DO JJ=1,G_NLOEN_MAX+R_NNOEXTZL+2
      DO JF=1,inf
         IGLG=D_NPTRLS(MYSETW)+KGL-1
         nlon = G_NLOEN(IGLG)
         !IOFF=D_NSTAGTF(KGL) ! modified
         if (JJ <= nlon) then
           IOFF=D_NSTAGTF(KGL)
           ZGTF(JF,IOFF+JJ)= ZGTF(JF,IOFF+JJ)/nlon
         end if

         IST  = 2*(G_NMEN(IGLG)+1)
         IF (ist+JJ <= nlon+R_NNOEXTZL+2) ZGTF(JF,IOFF+IST+JJ) = 0.0_JPRBT
         IF (nlon == 1) ZGTF(JF,IST+IOFF) = 0.0_JPRBT
      ENDDO
   ENDDO
ENDDO

write(nout,*) "transfer latitude buffer to FOUBUF_IN",kf_fs
! like: foubuf_in(1:2*kf_fs,ista) = zgtf(1:2*kf_fs,ioff+2*jm)

!$ACC PARALLEL LOOP DEFAULT(NONE) COLLAPSE(3) PRIVATE(IGLG,IPROC,ISTA,IOFF)
DO KGL=IBEG,IEND,IINC
   DO JM=0,G_NMEN_MAX
      DO JF=1,KF_FS
         IGLG = D_NPTRLS(MYSETW)+KGL-1
         if  (JM > G_NMEN(IGLG)) cycle

         IPROC = D_NPROCM(JM)
         ISTA  = (D_NSTAGT1B(D_MSTABF(IPROC))+D_NPNTGTB0(JM,KGL))*2*KF_FS
         IOFF  = 1+D_NSTAGTF(KGL)

         ! imaginary may be not JM+1 but JM+G_NMEN(IGLG)+1
         FOUBUF_IN(ISTA+2*JF-1) = ZGTF(2*JF-1,2*JM+IOFF)
         FOUBUF_IN(ISTA+2*JF) = ZGTF(2*JF,2*JM+IOFF)
      ENDDO
   ENDDO
END DO

!$ACC END DATA

#ifndef USE_CUDA_AWARE_MPI_FT
!$ACC UPDATE HOST(FOUBUF_IN)
#endif
END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
