      SUBROUTINE LIHEATT(CARTDT,DEPSVT,SGTOTT,DSTRAT,ELDIST,
     .                   GPCODT,LARGET,NDIMET,NDOFNT,NNODET,NSTR1T,
     .                   PROPST,SHAPET,STRANT,STRA0T,TSTRAT,XJACMT,
     .                   NDOFCT,SIGMAT,VELCMT,DMATXT,KDYNAT)
C***********************************************************************
C
C**** THIS ROUTINE EVALUATES TEMPERATURE GRADIENTS AND INCREMENTAL HEATS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'inte_omt.f'
C
      DIMENSION CARTDT(NDIMET,*), SGTOTT(*),        DSTRAT(*),
     .          ELDIST(NDOFCT,*), PROPST(*),
     .          GPCODT(*),        SHAPET(*),        STRANT(*),
     .          STRA0T(*),        TSTRAT(*),        XJACMT(NDIMET,*),
     .          SIGMAT(*),        VELCMT(*),        DMATXT(NSTR1T,*)
C
C**** CALCULATE THE TOTAL TEMPERATURE GRADIENT IN TSTRAT
C
      DO IDIMET=1,NDIMET
       TSTRAT(IDIMET)=0.0
       DO IEVABT=1,NNODET
        TSTRAT(IDIMET)=TSTRAT(IDIMET)+CARTDT(IDIMET,IEVABT)*
     .                 ELDIST(NDOFCT,IEVABT)
       END DO
      END DO
C
C**** CALCULATE THE INCREMENTAL TEMPERATURE GRADIENT IN DSTRAT
C
      DO IDIMET=1,NDIMET
       DSTRAT(IDIMET)=0.0
       IF(KDYNAT.EQ.1) THEN
        DO IEVABT=1,NNODET
         DSTRAT(IDIMET)=DSTRAT(IDIMET)+CARTDT(IDIMET,IEVABT)*
     .                  VELCMT(IEVABT)*DTIMET
        END DO
       ENDIF
      END DO
C
C**** CALCULATE THE EFFECTIVE HEATS 
C
      DO IDIMET=1,NDIMET
       SGTOTT(IDIMET)=0.0
       SIGMAT(IDIMET)=0.0
       DO JDIMET=1,NDIMET
        SGTOTT(IDIMET)=SGTOTT(IDIMET)+DMATXT(IDIMET,JDIMET)*
     .                 TSTRAT(JDIMET)
        SIGMAT(IDIMET)=SIGMAT(IDIMET)+DMATXT(IDIMET,JDIMET)*
     .                 DSTRAT(JDIMET)
       ENDDO
      ENDDO
C
      IF(KDYNAT.EQ.1) THEN         ! transient
       IF(KINTET.EQ.1) THEN        ! Euler's method
        DO IDIMET=1,NDIMET
         SGTOTT(IDIMET)=0.0
         SIGMAT(IDIMET)=0.0
         DO JDIMET=1,NDIMET
          SGTOTT(IDIMET)=SGTOTT(IDIMET)+DMATXT(IDIMET,JDIMET)*
     .                  (TALFAT*TSTRAT(JDIMET)+
     .                 (1.0D+00-TALFAT)*(TSTRAT(JDIMET)-DSTRAT(JDIMET)))
          SIGMAT(IDIMET)=SIGMAT(IDIMET)+DMATXT(IDIMET,JDIMET)*
     .                   DSTRAT(JDIMET)
         ENDDO
        ENDDO
       ENDIF
      ENDIF
C
C**** CALCULATES THE THERMAL DISIPATION
C
      auxtt=0.0
      auxti=0.0
      do idimet=1,ndimet
      auxtt=auxtt+sgtott(idimet)*tstrat(idimet)
      auxti=auxti+sigmat(idimet)*dstrat(idimet)
      end do
C
      dicai=auxti
      dicat=auxtt
C   
      RETURN
      END
