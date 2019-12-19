      SUBROUTINE SUMMARS
C***********************************************************************
C
C**** THIS ROUTINE WRITES A SUMMARY OF SPENT COMPUTING TIME
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'auxl_oms.f'
      INCLUDE 'prob_oms.f'
C
      DATA PERMXS,PERSTS,PERSFS,PERASS,PERSOS,PERRES,PERIOS,PERDAS,
     .     PEROUS,PERRSS,PERRNS,PEROTS/
     .         0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
     .         0.,    0.,    0.,    0./
C
      CALL CPUTIMS(CUCPUS)
C
      TOCPUS=CUCPUS-CPUINS
      CONTRS=CPURSS+CPUSTS+CPURES+CPUSFS+CPUASS+CPUSOS+CPUOUS+CPUDAS+
     .       CPURNS
      OTHERS=TOCPUS-CONTRS
      CPUMXS=CPUSTS+CPUSFS+CPUASS
      CPUIOS=CPUDAS+CPUOUS+CPURSS
C
      IF(TOCPUS.NE.0.0) PERMXS=100.0*CPUMXS/TOCPUS
      IF(CPUMXS.NE.0.0) THEN
        PERSTS=100.0*CPUSTS/CPUMXS
        PERSFS=100.0*CPUSFS/CPUMXS
        PERASS=100.0*CPUASS/CPUMXS
      ENDIF
      IF(TOCPUS.NE.0.0) THEN
       PERSOS=100.0*CPUSOS/TOCPUS
       PERRES=100.0*CPURES/TOCPUS
       PERIOS=100.0*CPUIOS/TOCPUS
      ENDIF
      IF(CPUIOS.NE.0.0) THEN
        PERDAS=100.0*CPUDAS/CPUIOS
        PEROUS=100.0*CPUOUS/CPUIOS
        PERRSS=100.0*CPURSS/CPUIOS
      ENDIF
      IF(TOCPUS.NE.0.0) THEN
       PERRNS=100.0*CPURNS/TOCPUS
       PEROTS=100.0*OTHERS/TOCPUS
      ENDIF
C
      WRITE(LURESS,900) TOCPUS,
     .               CPUMXS,PERMXS,CPUSTS,PERSTS,CPUSFS,PERSFS,CPUASS,
     .               PERASS,
     .               CPUSOS,PERSOS,CPURES,PERRES,
     .               CPUIOS,PERIOS,CPUDAS,PERDAS,CPUOUS,PEROUS,CPURSS,
     .               PERRSS,
     .               CPURNS,PERRNS,OTHERS,PEROTS
      IF(IEVFI.EQ.1)
     . WRITE(LUPRIS,900) TOCPUS,
     .               CPUMXS,PERMXS,CPUSTS,PERSTS,CPUSFS,PERSFS,CPUASS,
     .               PERASS,
     .               CPUSOS,PERSOS,CPURES,PERRES,
     .               CPUIOS,PERIOS,CPUDAS,PERDAS,CPUOUS,PEROUS,CPURSS,
     .               PERRSS,
     .               CPURNS,PERRNS,OTHERS,PEROTS
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
