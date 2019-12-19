      SUBROUTINE SUMMART
C***********************************************************************
C
C**** THIS ROUTINE WRITES A SUMMARY OF SPENT COMPUTING TIME
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_omt.f'
      INCLUDE 'prob_omt.f'
C
      DATA PERMXT,PERSTT,PERSFT,PERAST,PERSOT,PERRET,PERIOT,PERDAT,
     .     PEROUT,PERRST,PERRNT,PEROTT/
     .         0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
     .         0.,    0.,    0.,    0./
C
      CALL CPUTIMT(CUCPUT)
C
      TOCPUT=CUCPUT-CPUINT
      CONTRT=CPURST+CPUSTT+CPURET+CPUSFT+CPUAST+CPUSOT+CPUOUT+CPUDAT+
     .       CPURNT
      OTHERT=TOCPUT-CONTRT
      CPUMXT=CPUSTT+CPUSFT+CPUAST
      CPUIOT=CPUDAT+CPUOUT+CPURST
C
      IF(TOCPUT.NE.0.0) PERMXT=100.0*CPUMXT/TOCPUT
      IF(CPUMXT.NE.0.0) THEN
        PERSTT=100.0*CPUSTT/CPUMXT
        PERSFT=100.0*CPUSFT/CPUMXT
        PERAST=100.0*CPUAST/CPUMXT
      ENDIF
      IF(TOCPUT.NE.0.0) THEN
       PERSOT=100.0*CPUSOT/TOCPUT
       PERRET=100.0*CPURET/TOCPUT
       PERIOT=100.0*CPUIOT/TOCPUT
      ENDIF
      IF(CPUIOT.NE.0.0) THEN
        PERDAT=100.0*CPUDAT/CPUIOT
        PEROUT=100.0*CPUOUT/CPUIOT
        PERRST=100.0*CPURST/CPUIOT
      ENDIF
      IF(TOCPUT.NE.0.0) THEN
       PERRNT=100.0*CPURNT/TOCPUT
       PEROTT=100.0*OTHERT/TOCPUT
      ENDIF
C
      WRITE(LUREST,900) TOCPUT,
     .               CPUMXT,PERMXT,CPUSTT,PERSTT,CPUSFT,PERSFT,CPUAST,
     .               PERAST,
     .               CPUSOT,PERSOT,CPURET,PERRET,
     .               CPUIOT,PERIOT,CPUDAT,PERDAT,CPUOUT,PEROUT,CPURST,
     .               PERRST,
     .               CPURNT,PERRNT,OTHERT,PEROTT
      WRITE(LUPRIT,900) TOCPUT,
     .               CPUMXT,PERMXT,CPUSTT,PERSTT,CPUSFT,PERSFT,CPUAST,
     .               PERAST,
     .               CPUSOT,PERSOT,CPURET,PERRET,
     .               CPUIOT,PERIOT,CPUDAT,PERDAT,CPUOUT,PEROUT,CPURST,
     .               PERRST,
     .               CPURNT,PERRNT,OTHERT,PEROTT
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
