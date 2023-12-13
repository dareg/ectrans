! (C) Copyright 2000- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

INTERFACE
SUBROUTINE SPECNORM(PSPEC,KVSET,KMASTER,KRESOL,PMET,PNORM)

!**** *SPECNORM* - Compute global spectral norms

!     Purpose.
!     --------
!        Interface routine for computing spectral norms

!**   Interface.
!     ----------
!     CALL SPECNORM(...)

!     Explicit arguments : All arguments optional
!     --------------------
!     PSPEC(:,:)  - Spectral array
!     KVSET(:)    - "B-Set" for each field
!     KMASTER     - processor to recieve norms
!     KRESOL      - resolution tag  which is required ,default is the
!                   first defined resulution (input)
!     PMET(:)     - metric
!     PNORM(:)    - Norms (output for processor KMASTER)
!
!     Method.
!     -------

!     Externals.  SET_RESOL - set resolution
!     ----------  SPNORM_CTL - control routine

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

!ifndef INTERFACE

USE TPM_GEN         ,ONLY : NERR
!USE TPM_DIM
USE TPM_DISTR       ,ONLY : D, NPRTRV, MYSETV, MYPROC

USE SET_RESOL_MOD   ,ONLY : SET_RESOL
USE SPNORM_CTL_MOD  ,ONLY : SPNORM_CTL
USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS

!endif INTERFACE

IMPLICIT NONE

! Declaration of arguments


REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PSPEC(:,:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KVSET(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KMASTER
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(IN)  :: PMET(:)
REAL(KIND=JPRB)    ,OPTIONAL, INTENT(OUT) :: PNORM(:)
INTEGER(KIND=JPIM) ,OPTIONAL, INTENT(IN)  :: KRESOL
!ifndef INTERFACE
!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE SPECNORM

END INTERFACE
