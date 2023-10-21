!**** *LEDIR* - Direct Legendre transform.

! KM - zonal wavenumber
! KFC - number of field to transform
! PAIA - antisymmetric part of Fourier fields for zonal wavenumber KM
! PSIA - symmetric part of Fourier fields for zonal wavenumber KM
! POA1 -  spectral fields for zonal wavenumber KM

MODULE LEDIR_MOD
  USE PARKIND_ECTRANS  ,ONLY : JPIM ,JPIB    ,JPRB,  JPRBT
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
  USE TPM_DIM         ,ONLY : R, R_NDGNH,R_NSMAX,R_NTMAX
  USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
  USE TPM_FIELDS      ,ONLY : F, &
   	 & ZAA,ZAS,LDZAA,LDZAS,TDZAA,TDZAS,&
   	 & DZBST,DLDZBA,DLDZBS,DTDZBA,DTDZBS,&
   	 & DZCST,DZCAT,DLDZCA,DLDZCS,DTDZCA,DTDZCS,&
   	 & ZAMAX, ZSMAX,&
   	 & IF_FS_DIR,ZAA0,DZBST0,DZCAT0,ZAS0,DZCST0,KMLOC0
  USE TPM_DISTR
  USE TPM_GEN, ONLY: NOUT
  USE TPM_FLT
  USE BUTTERFLY_ALG_MOD
  USE CUDA_GEMM_BATCHED_MOD!!, ONLY: CUDA_TCGEMM_BATCHED, CUDA_GEMM_BATCHED
  USE CUBLAS_MOD, ONLY : CUDA_DGEMM_BATCHED
  USE, INTRINSIC :: ISO_C_BINDING
  USE IEEE_ARITHMETIC

  IMPLICIT NONE
CONTAINS
SUBROUTINE LEDIR2(KF_FS,KLED2,PAIA,POA1)
	INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS,KLED2
	REAL(KIND=JPRBT),    INTENT(IN)  :: PAIA(:,:,:)
	REAL(KIND=JPRBT),    INTENT(OUT) :: POA1(:,:,:)

	INTEGER(KIND=JPIM)  :: ISTAT,KM,KMLOC,KFC,KDGLU
	INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, IF, J, JK, IRET
	INTEGER(KIND=JPIM) :: ITHRESHOLD
	REAL(KIND=JPRB) :: RRPELTMDIR = 100.0_JPRB
	REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

	IF (LHOOK) CALL DR_HOOK('LE2_DGEMM',0,ZHOOK_HANDLE)

	KFC = 2*KF_FS

	!$ACC DATA PRESENT(F,F%RW,D,D_NUMP,D_MYMS,R,R_NDGNH,G,G_NDGLU,R_NSMAX,R_NTMAX,PAIA,&
	!$ACC ZAA,DZBST,DZCAT,POA1,dzbst0,dzcat0)

	! anti-symmetric

	!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISL,ISKIP) DEFAULT(NONE)
	DO KMLOC=1,D_NUMP
   	DO J=1,R_NDGNH   
      	DO JK=1,KFC
         	KM = D_MYMS(KMLOC)   
         	KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         	IF (J <= KDGLU) THEN
            	ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)

            	IF(KM == 0)THEN
               	ISKIP = 2
            	ELSE
               	ISKIP = 1
            	ENDIF

            	IF (MOD((JK-1),ISKIP) == 0) THEN
               	!DZBST((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*R_NDGNH)*IF_FS_DIR)=PAIA(JK,ISL+J-1,KMLOC)*F%RW(ISL+J-1)
               	DZBST((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZBA)*DTDZBA)=PAIA(JK,ISL+J-1,KMLOC)*F%RW(ISL+J-1)
            	END IF
         	END IF
      	ENDDO
   	ENDDO
	END DO

	! Get C in transpose format to get better memory access patterns later
	! C=A*B =>
	! C^T=B^T*A^T
	!$ACC HOST_DATA USE_DEVICE(ZAA,DZBST,DZCAT)
	CALL CUDA_GEMM_BATCHED('N', 'N',DTDZBA, TDZAA, DLDZBA,1.0_JPRBT,DZBST, DTDZBA, DLDZBA,&
	   ZAA, LDZAA, TDZAA,0._JPRBT,DZCAT, DTDZCA, DLDZCA,D_NUMP)
	!$ACC END HOST_DATA

	!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IA) DEFAULT(NONE)
	DO KMLOC=1,D_NUMP
   	DO J=1,(R_NTMAX+2)/2
      	DO JK=1,KFC
         	KM = D_MYMS(KMLOC)
         	IF(KM == 0)THEN
            	ISKIP = 2
         	ELSE
            	ISKIP = 1
         	ENDIF

         	IF (MOD((JK-1),ISKIP) == 0) THEN
            	ILA = (R_NTMAX-KM+2)/2
            	IA  = 1+MOD(R_NTMAX-KM+2,2)
            	IF (J <= ILA) THEN
               	POA1(JK,IA+(J-1)*2,KMLOC) = DZCAT((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZCA)*DTDZCA)
            	END IF
         	END IF
      	ENDDO
   	ENDDO
	ENDDO

	! compute m=0 in double precision:
	IF(KMLOC0 > 0) THEN
		write(nout,*) "ledir2 - kmloc0/ntmax/ndglu(0):",kmloc0,R_NTMAX,G_NDGLU(0)
   	ISKIP = 2
      KDGLU = MIN(R_NDGNH,G_NDGLU(0))
      ISL = MAX(R_NDGNH-G_NDGLU(0)+1,1)

	  !$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE)
	  DO J=1,KDGLU
   	 DO JK=1,KFC
         IF (MOD((JK-1),ISKIP) == 0) THEN
            DZBST0((JK-1)/ISKIP+1+(J-1)*DTDZBA)=PAIA(JK,ISL+J-1,KMLOC0)*F%RW(ISL+J-1)
         END IF
   	 ENDDO
	  ENDDO

	  !$ACC HOST_DATA USE_DEVICE(ZAA0,DZBST0,DZCAT0)
	  CALL CUDA_DGEMM_BATCHED('N','N',DTDZBA,int(TDZAA,kind=jpim),int(DLDZBA,kind=jpim), &
      	  & 1.0_JPRD,DZBST0,DTDZBA,int(DLDZBA,kind=jpim),&
      	  & ZAA0,LDZAA,int(TDZAA,kind=jpim),0._JPRD,DZCAT0,DTDZCA,int(DLDZCA,kind=jpim),1)
	  !$ACC END HOST_DATA

      ILA = (R_NTMAX+2)/2
      IA  = 1+MOD(R_NTMAX+2,2)

   	!$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE)
   	DO J=1,ila
   		DO JK=1,KFC
         	IF (MOD((JK-1),ISKIP) == 0) THEN
               POA1(JK,IA+(J-1)*2,KMLOC0) = DZCAT0((JK-1)/ISKIP+1+(J-1)*DTDZCA)
         	END IF
   		ENDDO
		ENDDO
	ENDIF

	!$ACC END DATA

	IF (LHOOK) CALL DR_HOOK('LE2_DGEMM',1,ZHOOK_HANDLE)
END SUBROUTINE

SUBROUTINE LEDIR1(KF_FS,KLED2,PAIA,POA1)
	INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS,KLED2
	REAL(KIND=JPRBT),    INTENT(IN)  :: PAIA(:,:,:)
	REAL(KIND=JPRBT),    INTENT(OUT) :: POA1(:,:,:)

	INTEGER(KIND=JPIM)  :: ISTAT,KM,KMLOC,KFC,KDGLU
	INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, IF, J, JK, IRET
	INTEGER(KIND=JPIM) :: ITHRESHOLD
	REAL(KIND=JPRB) :: RRPELTMDIR = 100.0_JPRB
	REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

	IF (LHOOK) CALL DR_HOOK('LE1_DGEMM',0,ZHOOK_HANDLE)

	KFC = 2*KF_FS

	write(nout,*) "ledir1 - nump/ndgnh/kfc:",D_NUMP,R_NDGNH,kfc

	!$ACC DATA PRESENT(F,F%RW,D,D_NUMP,D_MYMS,R,R_NDGNH,G,G_NDGLU,R_NSMAX,R_NTMAX,PAIA,&
	!$ACC ZAS,DZBST,DZCST,POA1,dzbst0,dzcst0)

	!$acc parallel loop collapse(3) private(KM,KDGLU,ISL,ISKIP) DEFAULT(NONE)
	DO KMLOC=1,D_NUMP
   	DO J=1,R_NDGNH   
      	DO JK=1,KFC
         	KM = D_MYMS(KMLOC)   
         	KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         	IF (J <= KDGLU) THEN
            	ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)

            	IF(KM == 0)THEN
               	ISKIP = 2
            	ELSE
               	ISKIP = 1
            	ENDIF

            	IF (MOD((JK-1),ISKIP) == 0) THEN
               	DZBST((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZBS)*DTDZBS)=PAIA(JK,ISL+J-1,KMLOC)*F%RW(ISL+J-1)
            	END IF
         	END IF
      	ENDDO
   	ENDDO
	END DO

	! Get C in transpose format to get better memory access patterns later
	! C=A*B =>
	! C^T=B^T*A^T
	!$ACC HOST_DATA USE_DEVICE(ZAS,DZBST,DZCST)
	CALL CUDA_GEMM_BATCHED('N', 'N',DTDZBS, TDZAS, DLDZBS,1.0_JPRBT,DZBST, DTDZBS, DLDZBS,&
	   ZAS, LDZAS, TDZAS,0._JPRBT,DZCST, DTDZCS, DLDZCS,D_NUMP)
	!$ACC END HOST_DATA

	!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IA,ILS,IS) DEFAULT(NONE)
	DO KMLOC=1,D_NUMP
   	DO J=1,(R_NTMAX+3)/2
      	DO JK=1,KFC
         	KM = D_MYMS(KMLOC)
         	IF(KM == 0)THEN
            	ISKIP = 2
         	ELSE
            	ISKIP = 1
         	ENDIF

         	IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
            	ILS = (R_NTMAX-KM+3)/2
            	IF (J <= ILS) THEN
               	IS  = 1+MOD(R_NTMAX-KM+1,2)
               	POA1(JK,IS+(J-1)*2,KMLOC) = DZCST((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZCS)*DTDZCS)            
            	END IF
         	END IF
      	ENDDO
   	ENDDO
	ENDDO

	IF(KMLOC0 > 0) THEN
		write(nout,*) "ledir1 - kmloc0/ntmax/ndglu(0):",kmloc0,R_NTMAX,G_NDGLU(0)
   	ISKIP = 2
      KDGLU = MIN(R_NDGNH,G_NDGLU(0))
      ISL = MAX(R_NDGNH-G_NDGLU(0)+1,1)

   	!$ACC PARALLEL LOOP COLLAPSE(2) DEFAULT(NONE)
   	DO J=1,KDGLU  
      	DO JK=1,KFC
         	IF (MOD((JK-1),ISKIP) == 0) THEN
               DZBST0((JK-1)/ISKIP+1+(J-1)*DTDZBS)=PAIA(JK,ISL+J-1,KMLOC0)*F%RW(ISL+J-1)
         	END IF
      	ENDDO
   	ENDDO

      !$ACC host_data use_device(ZAS0,DZBST0,DZCST0)
      call CUDA_DGEMM_BATCHED('N','N',DTDZBS,TDZAS,DLDZBS,1.0_JPRD,DZBST0,DTDZBS,DLDZBS,&
			ZAS0,LDZAS,TDZAS,0._JPRD,DZCST0,DTDZCS,DLDZCS,1)
      !$ACC end host_data

	   ILS = (R_NTMAX+3)/2
	   IS  = 1+MOD(R_NTMAX+1,2)

   	!$ACC parallel loop collapse(2) DEFAULT(NONE)
   	DO J=1,ils
      	DO JK=1,KFC
         	if (MOD(JK-1,ISKIP) == 0) then
            	POA1(JK,IS+(J-1)*2,KMLOC0) = DZCST0((JK-1)/ISKIP+1+(J-1)*DTDZCS)
         	end if
      	ENDDO
   	ENDDO
	ENDIF

	!$ACC END DATA

	IF (LHOOK) CALL DR_HOOK('LE1_DGEMM',1,ZHOOK_HANDLE)
END SUBROUTINE

SUBROUTINE LEDIR(KF_FS,KLED2,PAIA,POA1,KMODE)

INTEGER(KIND=JPIM), INTENT(IN)  :: KF_FS,KLED2
INTEGER(KIND=JPIM), INTENT(IN)  :: KMODE
REAL(KIND=JPRBT),    INTENT(IN)  :: PAIA(:,:,:)
REAL(KIND=JPRBT),    INTENT(OUT) :: POA1(:,:,:)

INTEGER(KIND=JPIM)  :: ISTAT,KM,KMLOC,KFC,KIFC,KDGLU
INTEGER(KIND=JPIM) :: IA, ILA, ILS, IS, ISKIP, ISL, IF, J, JK, IRET
INTEGER(KIND=JPIM) :: ITHRESHOLD
REAL(KIND=JPRB) :: RRPELTMDIR = 100.0_JPRB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',0,ZHOOK_HANDLE)

KFC = 2*KF_FS
KIFC = KFC

!$ACC DATA PRESENT(F,F%RW,D,D_NUMP,D_MYMS,R,R_NDGNH,G,G_NDGLU,R_NSMAX,R_NTMAX,PAIA,&
!$ACC ZAA,ZAS,DZBST,DZCST,DZCAT,POA1,dzbst0,dzcat0,dzcst0)

! anti-symmetric

IF ( KMODE == -1 ) THEN
	!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,KDGLU,ISL,ISKIP) DEFAULT(NONE)
	DO KMLOC=1,D_NUMP
   	DO J=1,R_NDGNH   
      	DO JK=1,KFC
         	KM = D_MYMS(KMLOC)   
         	KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         	IF (J .LE. KDGLU) THEN
            	ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)

            	IF(KM == 0)THEN
               	ISKIP = 2
            	ELSE
               	ISKIP = 1
            	ENDIF
            	IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
               	DZBST((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZBA)*DTDZBA)=PAIA(JK,ISL+J-1,KMLOC)*F%RW(ISL+J-1)
            	END IF
         	END IF
      	ENDDO
   	ENDDO
	END DO

	! Get C in transpose format to get better memory access patterns later
	!C=A*B =>
	! C^T=B^T*A^T
	!$ACC HOST_DATA USE_DEVICE(ZAA,DZBST,DZCAT)
	CALL CUDA_GEMM_BATCHED( &
	  & 'N', 'N', &
	  & DTDZBA, TDZAA, DLDZBA, &
	  & 1.0_JPRBT, &
	  & DZBST, DTDZBA, DLDZBA, &
	  & ZAA, LDZAA, TDZAA, &
	  & 0._JPRBT, &
	  & DZCAT, DTDZCA, DLDZCA, &
	  & D_NUMP)
	!$ACC END HOST_DATA

	!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IA,ILS) DEFAULT(NONE)
	DO KMLOC=1,D_NUMP
   	DO J=1,(R_NTMAX+2)/2
      	DO JK=1,KFC
         	KM = D_MYMS(KMLOC)
         	IF(KM == 0)THEN
            	ISKIP = 2
         	ELSE
            	ISKIP = 1
         	ENDIF

         	IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
            	ILA = (R_NTMAX-KM+2)/2
            	IA  = 1+MOD(R_NTMAX-KM+2,2)
            	IF (J .LE. ILA) THEN
               	POA1(JK,IA+(J-1)*2,KMLOC) = DZCAT((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZCA)*DTDZCA)
            	END IF
         	END IF
      	ENDDO
   	ENDDO
	ENDDO

	! compute m=0 in double precision:
	IF(KMLOC0 > 0) THEN
   	print*,'computing m=0 in double precision'
   	ISKIP = 2

	  !$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KM,KDGLU,ISL,ISKIP) DEFAULT(NONE)
	  DO J=1,R_NDGNH
   	 DO JK=1,KFC
         	KDGLU = MIN(R_NDGNH,G_NDGLU(0))
         	IF (J .LE. KDGLU) THEN
            	ISL = MAX(R_NDGNH-G_NDGLU(0)+1,1)
            	IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
               	DZBST0((JK-1)/ISKIP+1+(J-1)*DTDZBA)=PAIA(JK,ISL+J-1,KMLOC0)*F%RW(ISL+J-1)
            	END IF
         	END IF
   	 ENDDO
	  ENDDO


	  ! Get C in transpose format to get better memory access patterns later
	  !C=A*B =>
	  ! C^T=B^T*A^T

	  !$ACC HOST_DATA USE_DEVICE(ZAA0,DZBST0,DZCAT0)
	  CALL CUDA_DGEMM_BATCHED('N','N',DTDZBA,int(TDZAA,kind=jpim),int(DLDZBA,kind=jpim), &
      	  & 1.0_JPRD,DZBST0,DTDZBA,int(DLDZBA,kind=jpim),&
      	  & ZAA0,LDZAA,int(TDZAA,kind=jpim),0._JPRD,DZCAT0,DTDZCA,int(DLDZCA,kind=jpim),1)
	  !$ACC END HOST_DATA

   	!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(ILA,IA,ILS) DEFAULT(NONE)
   	DO J=1,(R_NTMAX+2)/2
   		DO JK=1,KFC
         	IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
            	ILA = (R_NTMAX+2)/2
            	IA  = 1+MOD(R_NTMAX+2,2)
            	IF (J .LE. ILA) THEN
               	POA1(JK,IA+(J-1)*2,KMLOC0) = DZCAT0((JK-1)/ISKIP+1+(J-1)*DTDZCA)
            	END IF
         	END IF
   		ENDDO
		ENDDO
	ENDIF
ELSE
	! symmetric

	!$acc parallel loop collapse(3) private(KM,KDGLU,ISL,ISKIP) DEFAULT(NONE)
	DO KMLOC=1,D_NUMP
   	DO J=1,R_NDGNH   
      	DO JK=1,KFC
         	KM = D_MYMS(KMLOC)   
         	KDGLU = MIN(R_NDGNH,G_NDGLU(KM))
         	IF (J .LE. KDGLU) THEN
            	ISL = MAX(R_NDGNH-G_NDGLU(KM)+1,1)

            	IF(KM == 0)THEN
               	ISKIP = 2
            	ELSE
               	ISKIP = 1
            	ENDIF
            	IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
	!               DZBST((JK-1)/ISKIP+1,J,KMLOC)=PSIA(JK,ISL+J-1,KMLOC)*F%RW(ISL+J-1)
               	DZBST((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZBS)*DTDZBS)=PAIA(JK,ISL+J-1,KMLOC)*F%RW(ISL+J-1)
            	END IF
         	END IF
      	ENDDO
   	ENDDO
	END DO

	! Get C in transpose format to get better memory access patterns later
	!C=A*B =>
	! C^T=B^T*A^T
	!$ACC HOST_DATA USE_DEVICE(ZAS,DZBST,DZCST)
	CALL CUDA_GEMM_BATCHED( &
	  & 'N', 'N', &
	  & DTDZBS, TDZAS, DLDZBS, &
	  & 1.0_JPRBT, &
	  & DZBST, DTDZBS, DLDZBS, &
	  & ZAS, LDZAS, TDZAS, &
	  & 0._JPRBT, &
	  & DZCST, DTDZCS, DLDZCS, &
	  & D_NUMP)
	!$ACC END HOST_DATA

	!$ACC PARALLEL LOOP COLLAPSE(3) PRIVATE(KM,ISKIP,ILA,IA,ILS,IS) DEFAULT(NONE)
	DO KMLOC=1,D_NUMP
   	DO J=1,(R_NTMAX+3)/2
      	DO JK=1,KFC

         	KM = D_MYMS(KMLOC)
         	IF(KM == 0)THEN
            	ISKIP = 2
         	ELSE
            	ISKIP = 1
         	ENDIF

         	IF (MOD((JK-1),ISKIP) .EQ. 0) THEN
            	ILS = (R_NTMAX-KM+3)/2
            	IF (J .LE. ILS) THEN
               	IS  = 1+MOD(R_NTMAX-KM+1,2)
               	POA1(JK,IS+(J-1)*2,KMLOC) = DZCST((JK-1)/ISKIP+1+(J-1+(KMLOC-1)*DLDZCS)*DTDZCS)            
            	END IF
         	END IF
      	ENDDO
   	ENDDO
	ENDDO

	IF(KMLOC0 > 0) THEN
	   ISKIP = 2
	   KDGLU = MIN(R_NDGNH,G_NDGLU(0))
      ISL = MAX(R_NDGNH-G_NDGLU(0)+1,1)

   	!$ACC PARALLEL LOOP COLLAPSE(2) PRIVATE(KDGLU,ISL) DEFAULT(NONE)
   	DO J=1,R_NDGNH   
      	DO JK=1,KFC
         	IF (J <= KDGLU) THEN
         	  IF (MOD((JK-1),ISKIP) == 0) THEN
               	DZBST0((JK-1)/ISKIP+1+(J-1)*DTDZBS)=PAIA(JK,ISL+J-1,KMLOC0)*F%RW(ISL+J-1)
         	  END IF
         	END IF
      	ENDDO
   	ENDDO

      !$ACC host_data use_device(ZAS0,DZBST0,DZCST0)
      call CUDA_DGEMM_BATCHED('N','N',DTDZBS,TDZAS,DLDZBS,1.0_JPRD,&
			DZBST0,DTDZBS,DLDZBS,ZAS0,LDZAS,TDZAS,0._JPRD,DZCST0,DTDZCS,DLDZCS,1)
      !$ACC end host_data

   	!$ACC parallel loop collapse(2) private(ILA,IA,ILS,IS) DEFAULT(NONE)
   	DO J=1,(R_NTMAX+3)/2
      	DO JK=1,KFC
         	if (MOD((JK-1),ISKIP) .eq. 0) then
            	ILS = (R_NTMAX+3)/2
            	if (J .le. ILS) then
               	IS  = 1+MOD(R_NTMAX+1,2)
               	POA1(JK,IS+(J-1)*2,KMLOC0) = DZCST0((JK-1)/ISKIP+1+(J-1)*DTDZCS)
            	end if
         	end if
      	ENDDO
   	ENDDO
	ENDIF
ENDIF

!$ACC END DATA

IF (LHOOK) CALL DR_HOOK('LE_DGEMM',1,ZHOOK_HANDLE)
END SUBROUTINE LEDIR
END MODULE LEDIR_MOD
