      SUBROUTINE MODUNS(WORKS, IP, MPRIN,
     .                  WORKST,IPT,MPRINT,
     .                  WORK1)
C=======================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8, ALLOCATABLE, INTENT(INOUT) :: WORK1(:)
C
C**** COUPLING VARIABLES
C
      INCLUDE 'nuec_om.f'
      INTERFACE
        INCLUDE 'modul0.inc'
        INCLUDE 'modul2.inc'
      END INTERFACE
C
C**** MECHANICAL VARIABLES (PERMANENT)
C
      DIMENSION WORKS(MPRIN)
      DIMENSION IP(50)                           ! see note in vul-*m*.f
C
C**** THERMAL VARIABLES (PERMANENT)
C
      DIMENSION WORKST(MPRINT), IPT(50)
C
C***********************************************************************
C
C**** ICOSTA=1 NON-CONVERGED STAGGERED STRATEGY
C
C***********************************************************************
C
C**** ITERME: DEFINES AN ASPECT OF THE STAGGERED STRATEGY USED
C
C              =0       THERMAL PROBLEM    >> MECHANICAL PROBLEM
C                       (CONSIDERING MECHANICAL INITIAL CONDITIONS -
C                        UNIDIRECTIONAL COUPLING)
C
C              =1       MECHANICAL PROBLEM >> THERMAL PROBLEM
C                       (WITHOUT CONSIDERING MECHANICAL INITIAL 
C                        CONDITIONS -
C                        BIDIRECTIONAL COUPLING)
C
C              =2       THERMAL PROBLEM    >> MECHANICAL PROBLEM
C                       (CONSIDERING MECHANICAL INITIAL CONDITIONS -
C                        BIDIRECTIONAL COUPLING)
C
C              =3       THERMAL PROBLEM      MECHANICAL PROBLEM
C                       (WITHOUT CONSIDERING MECHANICAL INITIAL 
C                        CONDITIONS -
C                        BIDIRECTIONAL COUPLING)
C
C              =4       THERMAL PROBLEM      MECHANICAL PROBLEM
C                       (CONSIDERING MECHANICAL INITIAL CONDITIONS -
C                        BIDIRECTIONAL COUPLING)
C
C              =5       THERMAL PROBLEM    >> MECHANICAL PROBLEM
C                       (WITHOUT CONSIDERING MECHANICAL INITIAL 
C                        CONDITIONS -
C                        BIDIRECTIONAL COUPLING)
C
C              =6       MECHANICAL PROBLEM >> THERMAL PROBLEM
C                       (CONSIDERING MECHANICAL INITIAL CONDITIONS -
C                        BIDIRECTIONAL COUPLING)
C
C***********************************************************************
C
C**** ISTAGG: DEFINES AN ASPECT STAGGERED STRATEGY USED
C
C     ISTAGG=0    >>    UNIDIRECTIONAL COUPLING IN THE CURRENT STEP
C                       (INCREMENTAL STAGGERED SOLUTION)
C
C     ISTAGG=1    >>    BIDIRECTIONAL COUPLING IN THE CURRENT STEP
C                       (ITERATIVE STAGGERED SOLUTION)
C
C***********************************************************************
C
      IF(ITERME.EQ.0) THEN
C
C**** UNIDIRECTIONAL COUPLING
C
       CALL RUNEND('moduns: such iterme not impl. yet')
       CALL MODUL0(
     .  WORKS(IP( 1)),WORKS(IP( 2)),WORKS(IP( 3)),WORKS(IP( 4)),
     .  WORKS(IP( 5)),WORKS(IP( 6)),WORKS(IP( 7)),WORKS(IP( 8)),
     .  WORKS(IP( 9)),WORKS(IP(10)),WORKS(IP(11)),WORKS(IP(12)),
     .  WORKS(IP(13)),WORKS(IP(14)),WORKS(IP(15)),WORKS(IP(16)),
     .  WORKS(IP(17)),WORKS(IP(18)),WORKS(IP(19)),WORKS(IP(20)),
     .  WORKS(IP(21)),WORKS(IP(22)),WORKS(IP(23)),WORKS(IP(24)),
     .  WORKS(IP(25)),WORKS(IP(26)),WORKS(IP(27)),WORKS(IP(28)),
     .  WORKS(IP(29)),WORKS(IP(30)),WORKS(IP(31)),WORKS(IP(32)),
     .  WORKS(IP(33)),WORKS(IP(34)),WORKS(IP(35)),WORKS(IP(36)),
     .  WORKS(IP(37)),
     .  WORK1,
     .  WORKST(IPT( 1)),WORKST(IPT( 2)),WORKST(IPT( 3)),WORKST(IPT( 4)),
     .  WORKST(IPT( 5)),WORKST(IPT( 6)),WORKST(IPT( 7)),WORKST(IPT( 8)),
     .  WORKST(IPT( 9)),WORKST(IPT(10)),WORKST(IPT(11)),WORKST(IPT(12)),
     .  WORKST(IPT(13)),WORKST(IPT(14)),WORKST(IPT(15)),WORKST(IPT(16)),
     .  WORKST(IPT(17)),WORKST(IPT(18)),WORKST(IPT(19)),WORKST(IPT(20)),
     .  WORKST(IPT(21)),WORKST(IPT(22)),WORKST(IPT(23)),WORKST(IPT(24)),
     .  WORKST(IPT(25)),WORKST(IPT(26)),WORKST(IPT(27)),WORKST(IPT(28)),
     .  WORKST(IPT(29)),WORKST(IPT(30)),WORKST(IPT(31)),
     .  WORK1)
C
       RETURN
C
      ENDIF
C
C**** BIDIRECTIONAL COUPLING
C
      GO TO (1,2,3,4,5,6), ITERME
C
    1 CONTINUE
      CALL RUNEND('moduns: such iterme not impl. yet')
      RETURN
C
    2 CONTINUE
      CALL RUNEND('moduns: such iterme not impl. yet')
      CALL MODUL2(
     .  WORKS(IP( 1)),WORKS(IP( 2)),WORKS(IP( 3)),WORKS(IP( 4)),
     .  WORKS(IP( 5)),WORKS(IP( 6)),WORKS(IP( 7)),WORKS(IP( 8)),
     .  WORKS(IP( 9)),WORKS(IP(10)),WORKS(IP(11)),WORKS(IP(12)),
     .  WORKS(IP(13)),WORKS(IP(14)),WORKS(IP(15)),WORKS(IP(16)),
     .  WORKS(IP(17)),WORKS(IP(18)),WORKS(IP(19)),WORKS(IP(20)),
     .  WORKS(IP(21)),WORKS(IP(22)),WORKS(IP(23)),WORKS(IP(24)),
     .  WORKS(IP(25)),WORKS(IP(26)),WORKS(IP(27)),WORKS(IP(28)),
     .  WORKS(IP(29)),WORKS(IP(30)),WORKS(IP(31)),WORKS(IP(32)),
     .  WORKS(IP(33)),WORKS(IP(34)),WORKS(IP(35)),WORKS(IP(36)),
     .  WORKS(IP(37)),
     .  WORK1,
     .  WORKST(IPT( 1)),WORKST(IPT( 2)),WORKST(IPT( 3)),WORKST(IPT( 4)),
     .  WORKST(IPT( 5)),WORKST(IPT( 6)),WORKST(IPT( 7)),WORKST(IPT( 8)),
     .  WORKST(IPT( 9)),WORKST(IPT(10)),WORKST(IPT(11)),WORKST(IPT(12)),
     .  WORKST(IPT(13)),WORKST(IPT(14)),WORKST(IPT(15)),WORKST(IPT(16)),
     .  WORKST(IPT(17)),WORKST(IPT(18)),WORKST(IPT(19)),WORKST(IPT(20)),
     .  WORKST(IPT(21)),WORKST(IPT(22)),WORKST(IPT(23)),WORKST(IPT(24)),
     .  WORKST(IPT(25)),WORKST(IPT(26)),WORKST(IPT(27)),WORKST(IPT(28)),
     .  WORKST(IPT(29)),WORKST(IPT(30)),WORKST(IPT(31)),
     .  WORK1)
      RETURN
C
    3 CONTINUE
      CALL RUNEND('moduns: such iterme not impl. yet')
      RETURN
C
    4 CONTINUE
      CALL RUNEND('moduns: such iterme not impl. yet')
      RETURN
C
    5 CONTINUE
      CALL RUNEND('moduns: such iterme not impl. yet')
      RETURN
C
    6 CONTINUE
      CALL RUNEND('moduns: such iterme not impl. yet')
      RETURN
C
      END
