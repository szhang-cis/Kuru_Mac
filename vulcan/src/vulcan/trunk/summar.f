      SUBROUTINE SUMMAR
C***********************************************************************
C
C**** THIS ROUTINE WRITES A SUMMARY OF SPENT COMPUTING TIME
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'
C
      DATA PERMX,PERST,PERSF,PERAS,PERSO,PERRE,PERIO,PERDA,
     .     PEROU,PERRS,PERRN,PEROT/
     .        0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
     .        0.,   0.,   0.,   0./
C
      CALL CPUTIM(CUCPU)
C
      TOCPU=CUCPU-CPUIN
      CONTR=CPURS+CPUST+CPURE+CPUSF+CPUAS+CPUSO+CPUOU+CPUDA+CPURN
      OTHER=TOCPU-CONTR
      CPUMX=CPUST+CPUSF+CPUAS
      CPUIO=CPUDA+CPUOU+CPURS
C
      IF(TOCPU.NE.0.0) PERMX=100.0*CPUMX/TOCPU
      IF(CPUMX.NE.0.0) THEN
        PERST=100.0*CPUST/CPUMX
        PERSF=100.0*CPUSF/CPUMX
        PERAS=100.0*CPUAS/CPUMX 
      ENDIF
      IF(TOCPU.NE.0.0) THEN
       PERSO=100.0*CPUSO/TOCPU
       PERRE=100.0*CPURE/TOCPU
       PERIO=100.0*CPUIO/TOCPU
      ENDIF
      IF(CPUIO.NE.0.0) THEN
        PERDA=100.0*CPUDA/CPUIO
        PEROU=100.0*CPUOU/CPUIO
        PERRS=100.0*CPURS/CPUIO
      ENDIF
      IF(TOCPU.NE.0.0) THEN
       PERRN=100.0*CPURN/TOCPU
       PEROT=100.0*OTHER/TOCPU
      ENDIF
C
      WRITE(LURES,900) TOCPU,
     .             CPUMX,PERMX,CPUST,PERST,CPUSF,PERSF,CPUAS,PERAS,
     .             CPUSO,PERSO,CPURE,PERRE,
     .             CPUIO,PERIO,CPUDA,PERDA,CPUOU,PEROU,CPURS,PERRS,
     .             CPURN,PERRN,OTHER,PEROT
      WRITE(LUPRI,900) TOCPU,
     .             CPUMX,PERMX,CPUST,PERST,CPUSF,PERSF,CPUAS,PERAS,
     .             CPUSO,PERSO,CPURE,PERRE,
     .             CPUIO,PERIO,CPUDA,PERDA,CPUOU,PEROU,CPURS,PERRS,
     .             CPURN,PERRN,OTHER,PEROT
C
      RETURN
  900 FORMAT(//,
     .  10X,'SUMMARY OF COMPUTING TIMES :',/,
     .  15X,'TOTAL CPU TIME             :',F10.2,/,
     .  15X,'COMPUTING MATRICES         :',F10.2,' (',F6.2,' % )',/,
     .  20X,'INITIAL COMPUTATIONS       :',F10.2,' (',F6.2,' % )',/,
     .  20X,'STIFFNESS MATRIX           :',F10.2,' (',F6.2,' % )',/,
     .  20X,'ASSEMBLY PROCESS           :',F10.2,' (',F6.2,' % )',/,
     .  15X,'SOLVING EQUATIONS          :',F10.2,' (',F6.2,' % )',/,
     .  15X,'EVALUATING RESIDUAL FORCES :',F10.2,' (',F6.2,' % )',/,
     .  15X,'I/O   OPERATIONS           :',F10.2,' (',F6.2,' % )',/,
     .  20X,'INPUT OPERATIONS           :',F10.2,' (',F6.2,' % )',/,
     .  20X,'OUTPUT TO EXTERNAL FILES   :',F10.2,' (',F6.2,' % )',/,
     .  20X,'OUTPUT TO RESTART  FILES   :',F10.2,' (',F6.2,' % )',/,
     .  15X,'RENUMBERING OPERATIONS     :',F10.2,' (',F6.2,' % )',/,
     .  15X,'OTHER OPERATIONS           :',F10.2,' (',F6.2,' % )')
      END
