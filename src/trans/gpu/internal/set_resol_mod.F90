! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SET_RESOL_MOD
CONTAINS
SUBROUTINE SET_RESOL(KRESOL,LDSETUP)
USE PARKIND1  ,ONLY : JPIM

USE TPM_GEN         ,ONLY : NOUT, MSETUP0, NCUR_RESOL, NMAX_RESOL,LENABLED
USE TPM_DIM         ,ONLY : R, DIM_RESOL,r_nsmax,r_ntmax,r_ndgl,r_ndgnh
!USE TPM_TRANS
USE TPM_DISTR       ,ONLY : D,DISTR_RESOL,d_nump,d_myms,d_nasm0,d_nptrls,d_nprocl,&
  d_nprocm,d_mstabf,d_nstagtf,d_nstagt0b,d_nstagt1b,d_npntgtb0,d_npntgtb1
USE TPM_GEOMETRY    ,ONLY : G,GEOM_RESOL,g_nloen,g_ndglu,g_nmen
USE TPM_FIELDS      ,ONLY : F, FIELDS_RESOL
USE TPM_FFT         ,ONLY : T, FFT_RESOL
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, FFTW_RESOL
#endif
USE TPM_FFTC        ,ONLY : TC, FFTC_RESOL
USE TPM_FLT
USE TPM_CTL        ,ONLY : C, CTL_RESOL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN) :: KRESOL
LOGICAL            ,OPTIONAL, INTENT(IN) :: LDSETUP

! Local varaibles
INTEGER(KIND=JPIM) :: IRESOL
LOGICAL :: LLSETUP

!     ------------------------------------------------------------------

IF (MSETUP0 == 0) CALL ABORT_TRANS('TRANS NOT SETUP')

LLSETUP = .FALSE.
IF(PRESENT(LDSETUP)) LLSETUP = LDSETUP

IRESOL = 1
IF (PRESENT(KRESOL)) THEN
  IF(KRESOL < 1 .OR. KRESOL > NMAX_RESOL) THEN
    CALL ABORT_TRANS('RESOL < 1 .OR. RESOL > NMAX_RESOL')
  ENDIF

  IRESOL = KRESOL
ENDIF

IF(.NOT.(LLSETUP.or.LENABLED(IRESOL))) CALL ABORT_TRANS('IRESOL NOT ENABLED OUT OF SETUP')

IF (IRESOL /= NCUR_RESOL) THEN
  write(nout,*) "--> changing resolution from/to:",ncur_resol,iresol
  NCUR_RESOL = IRESOL
  R => DIM_RESOL(NCUR_RESOL)
  F => FIELDS_RESOL(NCUR_RESOL)
  G => GEOM_RESOL(NCUR_RESOL)
  D => DISTR_RESOL(NCUR_RESOL)
  T => FFT_RESOL(NCUR_RESOL)
  TC => FFTC_RESOL(NCUR_RESOL)
  S => FLT_RESOL(NCUR_RESOL)
  C => CTL_RESOL(NCUR_RESOL)
#ifdef WITH_FFTW
  TW => FFTW_RESOL(NCUR_RESOL)
#endif

  if (.not.llsetup) then
	 r_nsmax = -1
	 !$acc update self(r_nsmax,r_ntmax,r_ndgl,r_ndgnh)
	 if (r%nsmax /= r_nsmax) call abort_trans("Error: r_nsmax")
	 if (r%ntmax /= r_ntmax) call abort_trans("Error: r_ntmax")
	 if (r%ndgl /= r_ndgl) call abort_trans("Error: r_ndgl")
	 if (r%ndgnh /= r_ndgnh) call abort_trans("Error: r_ndgnh")
	 !$acc update self(g_nloen,g_ndglu,g_nmen)
	 if (any(g%nloen /= g_nloen)) call abort_trans("Error: g_nloen")
	 if (any(g%ndglu /= g_ndglu)) call abort_trans("Error: g_ndglu")
	 if (any(g%nmen /= g_nmen)) call abort_trans("Error: g_nmen")

	 !$acc update self(d_nump,d_myms,d_nasm0,d_nptrls,d_nprocl,d_nprocm,d_mstabf,d_nstagtf,&
	 !$acc d_nstagt0b,d_nstagt1b,d_npntgtb0,d_npntgtb1)
	 if (d%nump /= d_nump) call abort_trans("Error: d_nump")
	 if (any(d%myms /= d_myms)) call abort_trans("Error: d_myms")
	 if (any(d%nasm0 /= d_nasm0)) call abort_trans("Error: d_nasm0")
	 if (any(d%nptrls /= d_nptrls)) call abort_trans("Error: d_nptrls")
	 if (any(d%nprocl /= d_nprocl)) call abort_trans("Error: d_nprocl")
	 if (any(d%nprocm /= d_nprocm)) call abort_trans("Error: d_nprocm")
	 if (any(d%mstabf /= d_mstabf)) call abort_trans("Error: d_mstabf")
	 if (any(d%nstagtf /= d_nstagtf)) call abort_trans("Error: d_nstagtf")
	 if (any(d%nstagt0b /= d_nstagt0b)) call abort_trans("Error: d_nstagt0b")
	 if (any(d%nstagt1b /= d_nstagt1b)) call abort_trans("Error: d_nstagt1b")
	 if (any(d%npntgtb0 /= d_npntgtb0)) call abort_trans("Error: d_npntgtb0")
	 if (any(d%npntgtb1 /= d_npntgtb1)) call abort_trans("Error: d_npntgtb1")
  end if
ENDIF

END SUBROUTINE SET_RESOL
END MODULE SET_RESOL_MOD
