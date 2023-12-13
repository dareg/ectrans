! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE TRANS_END(CDMODE)

!**** *TRANS_END* - Terminate transform package

!     Purpose.
!     --------
!     Terminate transform package. Release all allocated arrays.

!**   Interface.
!     ----------
!     CALL TRANS_END

!     Explicit arguments : None
!     --------------------

!     Method.
!     -------

!     Externals.  None
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!          G. Radnoti: 19-03-2009: intermediate end of transf to allow to switch to mono-task transforms
!        R. El Khatib 09-Jul-2013 LENABLED

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : MSETUP0, NCUR_RESOL, NMAX_RESOL, LENABLED,NDEF_RESOL
USE TPM_DIM         ,ONLY : R, DIM_RESOL
USE TPM_DISTR       ,ONLY : D, DISTR_RESOL, NPRCIDS
USE TPM_GEOMETRY    ,ONLY : G, GEOM_RESOL
USE TPM_FIELDS      ,ONLY : F, FIELDS_RESOL
USE TPM_FFT         ,ONLY : T, FFT_RESOL, TB, FFTB_RESOL
USE TPM_CTL         ,ONLY : C, CTL_RESOL
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, FFTW_RESOL
#endif
USE TPM_FLT
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN

USE EQ_REGIONS_MOD  ,ONLY : N_REGIONS
USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE DEALLOC_RESOL_MOD   ,ONLY : DEALLOC_RESOL
!

IMPLICIT NONE
CHARACTER*5, OPTIONAL,  INTENT(IN) :: CDMODE
! Local variables
!endif INTERFACE

END SUBROUTINE TRANS_END
END INTERFACE
