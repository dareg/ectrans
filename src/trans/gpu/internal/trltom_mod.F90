! (C) Copyright 1995- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE TRLTOM_MOD
  CONTAINS
  SUBROUTINE TRLTOM_CUDAAWARE(PFBUF_IN,PFBUF,KFIELD)
  
  !**** *TRLTOM * - transposition in Fourierspace
  
  !     Purpose.
  !     --------
  !              Transpose Fourier coefficients from partitioning
  !              over latitudes to partitioning over wave numbers
  !              This is done between inverse Legendre Transform
  !              and inverse FFT.
  !              This is the inverse routine of TRMTOL.
  
  !**   Interface.
  !     ----------
  !        *CALL* *TRLTOM(...)*
  
  !        Explicit arguments : PFBUF  - Fourier coefficient buffer. It is
  !        --------------------          used for both input and output.
  
  !                             KFIELD - Number of fields communicated
  
  !        Implicit arguments :
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.
  !     ----------
  
  !     Reference.
  !     ----------
  !        ECMWF Research Department documentation of the IFS
  
  !     Author.
  !     -------
  !        MPP Group *ECMWF*
  
  !     Modifications.
  !     --------------
  !        Original : 95-10-01
  !        Modified : 97-06-18 G. Mozdzynski - control MPI mailbox use
  !                                            (NCOMBFLEN) for nphase.eq.1
  !        Modified : 99-05-28  D.Salmond - Optimise copies.
  !        Modified : 00-02-02  M.Hamrud  - Remove NPHASE
  !        D.Salmond : 01-11-23 LIMP_NOOLAP Option for non-overlapping message
  !                             passing and buffer packing
  !        G.Mozdzynski : 08-01-01 Cleanup
  !        Y.Seity   : 07-08-30 Add barrier synchonisation under LSYNC_TRANS
  !     ------------------------------------------------------------------
  
  USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
  
  USE MPL_MODULE  ,ONLY : MPL_ALLTOALLV, MPL_BARRIER, MPL_ALL_MS_COMM, MPL_MYRANK, MPL_WAIT, JP_NON_BLOCKING_STANDARD
  
  USE TPM_DISTR       ,ONLY : D, MTAGLM, MYSETW, NPRTRW, NPROC, MYPROC
  USE TPM_GEN         ,ONLY : LSYNC_TRANS
  
  USE MPI
  
  !USE SET2PE_MOD
  !USE MYSENDSET_MOD
  !USE MYRECVSET_MOD
  !USE ABORT_TRANS_MOD
  !
  
  IMPLICIT NONE
  
  
  INTERFACE
  
    FUNCTION ALLTOALLV_CUDAIPC(input,len,soff,output,roff,mtol_or_ltom) BIND(C,name='Alltoallv_CUDAIPC')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      real(c_double), dimension(*) :: input,output
      integer(c_int), dimension(*) :: len,soff,roff
      integer(c_int),value :: mtol_or_ltom
      integer(c_int) :: ALLTOALLV_CUDAIPC
    END FUNCTION ALLTOALLV_CUDAIPC
  
  END INTERFACE
  
  
  INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
  REAL(KIND=JPRBT)   ,INTENT(INOUT)  :: PFBUF(:)
  REAL(KIND=JPRBT)   ,INTENT(INOUT)  :: PFBUF_IN(:)
  
  INTEGER(KIND=JPIM) :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)
  
  INTEGER(KIND=JPIM) :: ITAG, J, ILEN, ISTA
  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR2
  
  REAL(KIND=JPRBT)    :: ZDUM(1)
  INTEGER(KIND=JPIM) :: IREQ
  INTEGER(KIND=JPIM) :: IERROR
  !     ------------------------------------------------------------------
  
  REAL(KIND=JPRBT) :: T1, T2, TIMEF, Tc
  INTEGER(KIND=JPIM) :: MTOL_OR_LTOM, NOFULLPEERACCESS
  INTEGER(KIND=JPIM) :: IRANK,IUNIT
  INTEGER(KIND=JPIM) :: FROM_SEND,FROM_RECV,TO_RECV,TO_SEND


  IF (LHOOK) CALL DR_HOOK('TRLTOM_CUDAAWARE',0,ZHOOK_HANDLE)

#ifdef PARKINDTRANS_SINGLE
#define TRLTOM_DTYPE MPI_REAL
#else
#define TRLTOM_DTYPE MPI_DOUBLE_PRECISION
#endif
  
  ITAG = MTAGLM
  
  DO J=1,NPRTRW
    ILENS(J) = D%NLTSGTB(J)*KFIELD
    IOFFS(J) = D%NSTAGT1B(D%MSTABF(J))*KFIELD
    ILENR(J) = D%NLTSFTB(J)*KFIELD
    IOFFR(J) = D%NSTAGT1B(J)*KFIELD
  ENDDO
  
  IF(NPROC > 1) THEN
    IF (LSYNC_TRANS) THEN
      CALL MPL_BARRIER(CDSTRING='TRLTOM BARRIER')
    ENDIF

    ! copy to self workaround
    IRANK = MPL_MYRANK(MPL_ALL_MS_COMM)
    IF (ILENS(IRANK) .ne. ILENR(IRANK)) THEN
        PRINT *, "ERROR", ILENS(IRANK), ILENR(IRANK)
        stop 1
    ENDIF

    IF (ILENS(IRANK) > 0) THEN
        FROM_SEND = IOFFS(IRANK) + 1
        TO_SEND = FROM_SEND + ILENS(IRANK) - 1
        FROM_RECV = IOFFR(IRANK) + 1
        TO_RECV = FROM_RECV + ILENR(IRANK) - 1
        !$ACC KERNELS ASYNC(1) DEFAULT(NONE) PRESENT(PFBUF,PFBUF_IN)
        PFBUF(FROM_RECV:TO_RECV) = PFBUF_IN(FROM_SEND:TO_SEND)
        !$ACC END KERNELS
        ILENS(IRANK) = 0
        ILENR(IRANK) = 0
    ENDIF
  
    !$ACC HOST_DATA USE_DEVICE(PFBUF_IN,PFBUF)

    CALL MPI_ALLTOALLV(PFBUF_IN,ILENS,IOFFS,TRLTOM_DTYPE,PFBUF,ILENR,IOFFR,TRLTOM_DTYPE,&
      MPL_ALL_MS_COMM,IERROR)

    !$ACC END HOST_DATA
    !$ACC WAIT(1)

    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  ELSE
    ILEN = D%NLTSGTB(MYSETW)*KFIELD
    ISTA = D%NSTAGT1B(MYSETW)*KFIELD

    !$ACC PARALLEL LOOP DEFAULT(NONE) PRESENT(PFBUF,PFBUF_IN)
    DO J=ISTA+1,ISTA+ILEN
      PFBUF(J) = PFBUF_IN(J)
    ENDDO
  ENDIF
  
  IF (LHOOK) CALL DR_HOOK('TRLTOM_CUDAAWARE',1,ZHOOK_HANDLE)
  !     ------------------------------------------------------------------
  END SUBROUTINE TRLTOM_CUDAAWARE

  SUBROUTINE TRLTOM(PFBUF_IN,PFBUF,KFIELD)
  
  !**** *TRLTOM * - transposition in Fourierspace
  
  !     Purpose.
  !     --------
  !              Transpose Fourier coefficients from partitioning
  !              over latitudes to partitioning over wave numbers
  !              This is done between inverse Legendre Transform
  !              and inverse FFT.
  !              This is the inverse routine of TRMTOL.
  
  !**   Interface.
  !     ----------
  !        *CALL* *TRLTOM(...)*
  
  !        Explicit arguments : PFBUF  - Fourier coefficient buffer. It is
  !        --------------------          used for both input and output.
  
  !                             KFIELD - Number of fields communicated
  
  !        Implicit arguments :
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.
  !     ----------
  
  !     Reference.
  !     ----------
  !        ECMWF Research Department documentation of the IFS
  
  !     Author.
  !     -------
  !        MPP Group *ECMWF*
  
  !     Modifications.
  !     --------------
  !        Original : 95-10-01
  !        Modified : 97-06-18 G. Mozdzynski - control MPI mailbox use
  !                                            (NCOMBFLEN) for nphase.eq.1
  !        Modified : 99-05-28  D.Salmond - Optimise copies.
  !        Modified : 00-02-02  M.Hamrud  - Remove NPHASE
  !        D.Salmond : 01-11-23 LIMP_NOOLAP Option for non-overlapping message
  !                             passing and buffer packing
  !        G.Mozdzynski : 08-01-01 Cleanup
  !        Y.Seity   : 07-08-30 Add barrier synchonisation under LSYNC_TRANS
  !     ------------------------------------------------------------------
  
  USE PARKIND_ECTRANS  ,ONLY : JPIM     ,JPRBT
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
  
  USE MPL_MODULE  ,ONLY : MPL_ALLTOALLV, MPL_BARRIER, MPL_ALL_MS_COMM, MPL_MYRANK, MPL_WAIT, JP_NON_BLOCKING_STANDARD
  
  USE TPM_DISTR       ,ONLY : D, MTAGLM, MYSETW, NPRTRW, NPROC, MYPROC
  USE TPM_GEN         ,ONLY : LSYNC_TRANS
  
  USE MPI
  
  !USE SET2PE_MOD
  !USE MYSENDSET_MOD
  !USE MYRECVSET_MOD
  !USE ABORT_TRANS_MOD
  !
  
  IMPLICIT NONE
  
  
  INTERFACE
  
    FUNCTION ALLTOALLV_CUDAIPC(input,len,soff,output,roff,mtol_or_ltom) BIND(C,name='Alltoallv_CUDAIPC')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      real(c_double), dimension(*) :: input,output
      integer(c_int), dimension(*) :: len,soff,roff
      integer(c_int),value :: mtol_or_ltom
      integer(c_int) :: ALLTOALLV_CUDAIPC
    END FUNCTION ALLTOALLV_CUDAIPC
  
  END INTERFACE
  
  
  INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
  REAL(KIND=JPRBT)   ,INTENT(INOUT)  :: PFBUF(:)
  REAL(KIND=JPRBT)   ,INTENT(INOUT)  :: PFBUF_IN(:)
  
  INTEGER(KIND=JPIM) :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)
  
  INTEGER(KIND=JPIM) :: ITAG, J, ILEN, ISTA
  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE_BAR2
  
  REAL(KIND=JPRBT)    :: ZDUM(1)
  INTEGER(KIND=JPIM) :: IREQ
  INTEGER(KIND=JPIM) :: IERROR
  !     ------------------------------------------------------------------
  
  REAL(KIND=JPRBT) :: T1, T2, TIMEF, tc
  INTEGER(KIND=JPIM) :: MTOL_OR_LTOM, NOFULLPEERACCESS
  INTEGER(KIND=JPIM) :: IRANK,iunit

  IF (LHOOK) CALL DR_HOOK('TRLTOM',0,ZHOOK_HANDLE)
  
  ITAG = MTAGLM
  
  DO J=1,NPRTRW
    ILENS(J) = D%NLTSGTB(J)*KFIELD
    IOFFS(J) = D%NSTAGT1B(D%MSTABF(J))*KFIELD
    ILENR(J) = D%NLTSFTB(J)*KFIELD
    IOFFR(J) = D%NSTAGT1B(J)*KFIELD
  ENDDO
  
  IF(NPROC > 1) THEN
    CALL MPL_ALLTOALLV(PSENDBUF=PFBUF_IN,KSENDCOUNTS=ILENS,&
     & PRECVBUF=PFBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
     & KCOMM=MPL_ALL_MS_COMM,CDSTRING='TRLTOM:')

    call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  ELSE
    ILEN = D%NLTSGTB(MYSETW)*KFIELD
    ISTA = D%NSTAGT1B(MYSETW)*KFIELD

    DO J=ISTA+1,ISTA+ILEN
      PFBUF(J) = PFBUF_IN(J)
    ENDDO
  ENDIF
  
  IF (LHOOK) CALL DR_HOOK('TRLTOM',1,ZHOOK_HANDLE)
  !     ------------------------------------------------------------------
  END SUBROUTINE TRLTOM
  END MODULE TRLTOM_MOD
