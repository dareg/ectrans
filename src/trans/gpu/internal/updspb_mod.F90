MODULE UPDSPB_MOD
CONTAINS
SUBROUTINE UPDSPB(KFIELD,POA,PSPEC,KFLDPTR)
  USE PARKIND_ECTRANS,ONLY : JPIM,JPRB,JPRBT
  USE TPM_DIM,ONLY : R,R_NSMAX,R_NTMAX
  USE TPM_DISTR,ONLY : D,D_NUMP,D_MYMS,D_NASM0

  IMPLICIT NONE

  INTEGER(KIND=JPIM),INTENT(IN) :: KFIELD
  REAL(KIND=JPRBT),INTENT(IN) :: POA(:,:,:)
  REAL(KIND=JPRB),INTENT(OUT) :: PSPEC(:,:)
  INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KFLDPTR(:)

  INTEGER(KIND=JPIM) :: KM,KMLOC,INM,JFLD,JN,IASM0

  ! The following transfer reads :
  ! SPEC(k,NASM0(m)+NLTN(n)*2)  =POA(nn,2*k-1) (real part)
  ! SPEC(k,NASM0(m)+NLTN(n)*2+1)=POA(nn,2*k  ) (imaginary part)
  ! with n from m to NSMAX
  ! and nn=NTMAX+2-n from NTMAX+2-m to NTMAX+2-NSMAX.
  ! NLTN(m)=NTMAX+2-m : n=NLTN(nn),nn=NLTN(n)
  ! nn is the loop index.

  IF (PRESENT(KFLDPTR)) stop 'Error: code path not (yet) supported in GPU version'

  !$ACC DATA PRESENT(PSPEC,POA,R,D)

  !$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,IASM0,INM) DEFAULT(NONE)
  DO KMLOC=1,D%NUMP     
     DO JN=0,R%NSMAX
        DO JFLD=1,KFIELD
           KM = D%MYMS(KMLOC)
           IASM0 = D%NASM0(KM)

           IF (KM == 0) THEN
              INM = IASM0+2*JN
              PSPEC(JFLD,INM) = POA(2*JFLD-1,R%NTMAX+2-JN,KMLOC)
              PSPEC(JFLD,INM+1) = 0.0_JPRBT
           ELSE
              if (JN >= KM) then
                 INM = IASM0+2*(JN-KM)
                 PSPEC(JFLD,INM) = POA(2*JFLD-1,R%NTMAX+2-JN,KMLOC)
                 PSPEC(JFLD,INM+1) = POA(2*JFLD,R%NTMAX+2-JN,KMLOC)
              end if
           end if
        ENDDO
     ENDDO
  ENDDO
  !$ACC END PARALLEL

  !$ACC END DATA
END SUBROUTINE UPDSPB
END MODULE UPDSPB_MOD
