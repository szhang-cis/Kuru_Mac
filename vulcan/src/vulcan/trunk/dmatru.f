      SUBROUTINE DMATRU(EPMTX,NSTRS,NTYPE,PROPS)
C************************************************************************
C
C****THIS ROUTINE SETS UP THE ELASTIC MATRIX FOR PLANE STRESS & STRAIN
C    AXISYMMETRIC AND 3-D joints elements
C
C************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION EPMTX(*), PROPS(*)
C
C***INITIALISE THE ELASTIC MATRIX
C
      NKOST=(NSTRS+1)*NSTRS/2
      DO 10 IKOST=1,NKOST
   10 EPMTX(IKOST)=0.0D00
C
C***READ ELASTIC PROPERTIES
C
      EXX=PROPS(2)
      PXY=PROPS(3)
      EYY=PROPS(4)
      PXZ=PROPS(5)
      EZZ=PROPS(6)
      PYZ=PROPS(7)
      GXY=PROPS(8)
      GXZ=PROPS(9)
      GYZ=PROPS(10)
C
      GOTO (1,2,3,4,5),NTYPE
C
C***PLANE STRESS
C
    1 EPMTX(1) =EXX
      EPMTX(2) =0.0
      EPMTX(4) =EYY
      EPMTX(6) =GXY
C
      RETURN
C
C***PLANE STRAIN
C
    2 EPMTX(1) =EXX
      EPMTX(2) =0.0
      EPMTX(4) =EXX*PXY
      EPMTX(5) =EYY
      EPMTX(7) =0.0
      EPMTX(8) =GXY
      EPMTX(10)=EZZ
C
      RETURN
C
C***AXI-SYMMETRIC CASE
C
    3 EPMTX(1) =EXX
      EPMTX(2) =0.0
      EPMTX(4) =EXX*PXY
      EPMTX(5) =EYY
      EPMTX(7) =0.0
      EPMTX(8) =GXY
      EPMTX(10)=EZZ
C
      RETURN
C
C***3-D CASE
C
c
c    not implemented yet!!!!!!!!!   
c
    4 CO1=1-PYZ*PYZ*EYY/EZZ
      CO2=PXY+PXZ*PYZ*EYY/EZZ
      C21=CO2/CO1
      CO3=1-(PXZ*PXZ+(EZZ/EYY)*C21*CO2)*EXX/EZZ
      CO4=PXZ+PYZ*C21
C
      EPMTX(1) =EXX/CO3
      EPMTX(2) =EPMTX(1)*C21
      EPMTX(4) =EPMTX(1)*CO4
      EPMTX(7) =EPMTX(2)*C21+EYY/CO1
      EPMTX(9) =EPMTX(4)*C21+EYY*PYZ/CO1
      EPMTX(12)=GXY
      EPMTX(16)=EPMTX(4)*PXZ+EPMTX(9)*PYZ+EZZ
      EPMTX(19)=GXZ
      EPMTX(21)=GYZ
C
      RETURN
C
C***1-D CASE
C
    5 EPMTX(1) =EXX
C
      RETURN
      END     
