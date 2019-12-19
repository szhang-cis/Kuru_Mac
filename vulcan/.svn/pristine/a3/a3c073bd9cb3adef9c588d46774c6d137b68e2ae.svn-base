      SUBROUTINE PROPANI(ITAPE,PROPS,VANIS,MATNO,PROEL,COORD,LNODS)
C***********************************************************************
C
C**** THIS ROUTINE READS THE MATERIAL SYSTEM OF COORDINATES FOR
C     ANISOTROPIC MODELS
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C**** MECHANICAL VARIABLES
C
      INCLUDE 'auxl_om.f'
      INCLUDE 'prob_om.f'
      INCLUDE 'inte_om.f'
      INCLUDE 'inpo_om.f'
C
      DIMENSION PROPS(NPROP,*),     VANIS(NANIV,NANIC,NELEM),
     .          MATNO(NELEM),       PROEL(NPREL,NGRUP),
     .          COORD(NDIME,NPOIN), LNODS(NNODE,NELEM)
      DIMENSION VECTD(3), COORO(3), COORP(3), DCRUP(3)
C
C**** READS 'MATERIAL_SYSTEM_OF_COORDINATES' FOR ANISOTROPY
C
C     Note for IFREN=2,3:
C           For planar anisotropy (normal isotropy), only one vector
C           must be input in 2D while two vectors must be input in 3D.
C           2D: v1 (input), v3 (out-of-plane) => v2=v3 x v1
C           3D: v1 (input), v2 (input)        => v3=v1 x v2
C           For full anisotropy, two vectors must be input in 3D.
C           3D: v1 (input), v2 (input)        => v3=v1 x v2
C
C     Notes for IFREN=55,57:
C           Input axis direction, axis point and angle
C           The results are independent of the coordinates of the axis
C           point
C
C     Note for IFREN=56:
C           Input fiber and surface (3D) vectors
C
      IF(IANIS.EQ.1) THEN
       NPRIN=0
       CALL LISTEN('PROPANI',NPRIN,ITAPE)
       IF(WORDS(1).EQ.'MATER') THEN
        IVECT=0                            ! vector d form
        IF(WORDS(2).EQ.'VECTO') IVECT=1    ! vector form
        CALL LISTEN('PROPANI',NPRIN,ITAPE)
        NELEA=INT(PARAM(1))
        IF(NELEA.LE.0.OR.NELEA.GT.NELEM)
     .   CALL RUNEND('PROPANI: WRONG NUMBER OF ANISOTROPIC ELEMENTS')
C
        DO IELEM=1,NELEM                   ! initialization
         DO IDIM1=1,NANIV
          DO IDIM2=1,NANIC
           VANIS(IDIM1,IDIM2,IELEM)=0.0D0
          END DO
         END DO
        END DO
C
        DO IELEA=1,NELEA
         CALL LISTEN('PROPANI',NPRIN,ITAPE)
         IELEX=INT(PARAM(1))
         LGRUP=MATNO(IELEX)                ! set number
         LMATS=INT(PROEL(1,LGRUP))         ! material number
         NNODL=INT(PROEL(2,LGRUP))
         LTYPE=INT(PROEL(5,LGRUP))
         IF(LTYPE.EQ.32)
     .    CALL RUNEND('PROPANI: CONTACT ELEMENT IN MATER._SYST._COORD.')
         CALL IDENPR(PROPS(1,LMATS))       ! obtains material properties
         ISOTR=INT(PROPS(1,LMATS))         ! 0=isotropic; 1=orthotropic
         IF(ISOTR.EQ.1) THEN
          IF(NFIBN.EQ.0) THEN
           IF(IFREN.NE.2.AND.IFREN.NE.3.AND.
     .        IFREN.NE.55.AND.IFREN.NE.56.AND.IFREN.NE.57)
     .      CALL RUNEND('PROPANI: ONLY IFREN=2,3,55-57 ARE AVAILABLE')
          ELSE
           IF(IFREN.NE.2.AND.IFREN.NE.3.AND.
     .        IFREN.NE.51.AND.IFREN.NE.52.AND.IFREN.NE.53.AND.
     .        IFREN.EQ.54.AND.IFREN.NE.58)
     .      CALL RUNEND('PROPANI:ONLY IFREN=2,3,51-54,58 ARE AVAILABLE')
          ENDIF
          IF(IFREN.EQ.55.OR.               ! Holzapfel model
     .       IFREN.EQ.57) THEN             ! Gasser model
           IF(IVECT.EQ.0) THEN             ! vector d form
            IF(NTYPE.EQ.3) THEN            ! axisymmetry
             VECTD(1)=PARAM(2)             ! vector d (// to axis)
             VECTD(2)=PARAM(3)
             VECTD(3)=0.0D0
             IF(PARAM(2).EQ.0.0D0.AND.PARAM(3).EQ.0.0D0) VECTD(3)=1.0D0
             COORO(1)=PARAM(4)             ! coordinates of axis point
             COORO(2)=PARAM(5)
             COORO(3)=0.0D0
             ANGLD=PARAM(6)                ! angle phi
            ENDIF
            IF(NTYPE.EQ.4) THEN            ! 3D
             VECTD(1)=PARAM(2)             ! vector d (// to axis)
             VECTD(2)=PARAM(3)
             VECTD(3)=PARAM(4)
             COORO(1)=PARAM(5)             ! coordinates of axis point
             COORO(2)=PARAM(6)
             COORO(3)=PARAM(7)
             ANGLD=PARAM(8)                ! angle phi
            ENDIF
C
            DO IDIME=1,3                   ! computation of point P
             COORP(IDIME)=0.0D0
            ENDDO
            DO INODL=1,NNODL
             IPOIN=LNODS(INODL,IELEX)
             DO IDIME=1,NDIME
              COORP(IDIME)=COORP(IDIME)+COORD(IDIME,IPOIN)
             ENDDO
            ENDDO
            DO IDIME=1,NDIME               ! average of nodal coord.
             COORP(IDIME)=COORP(IDIME)/NNODL
            ENDDO
C
            DCRUP(1)=VECTD(2)*(COORP(3)-COORO(3))-   ! d x OP
     .               VECTD(3)*(COORP(2)-COORO(2))
            DCRUP(2)=VECTD(3)*(COORP(1)-COORO(1))-
     .               VECTD(1)*(COORP(3)-COORO(3))
            DCRUP(3)=VECTD(1)*(COORP(2)-COORO(2))-
     .               VECTD(2)*(COORP(1)-COORO(1))
C
            DCRUM=DCRUP(1)*DCRUP(1)+DCRUP(2)*DCRUP(2)+DCRUP(3)*DCRUP(3)
            DCRUM=DSQRT(DCRUM)                       ! modulus of d x OP
            VECTM=VECTD(1)*VECTD(1)+VECTD(2)*VECTD(2)+VECTD(3)*VECTD(3)
            VECTM=DSQRT(VECTM)                       ! modulus of d
C
            TWOPI=6.283185307179586D0
            ANGLD=ANGLD*TWOPI/360.0D0
C
            DO IDIME=1,3
             VANIS(1,IDIME,IELEX)=VECTD(IDIME)/VECTM*DCOS(ANGLD)+  ! a0
     .                            DCRUP(IDIME)/DCRUM*DSIN(ANGLD)
             VANIS(2,IDIME,IELEX)=VECTD(IDIME)/VECTM*DCOS(ANGLD)-  ! b0
     .                            DCRUP(IDIME)/DCRUM*DSIN(ANGLD)
            ENDDO
           ENDIF
C
           IF(IVECT.EQ.1) THEN                       ! vector form
            DO IDIME=1,3
             VANIS(1,IDIME,IELEX)=PARAM(1+IDIME)                   ! a0
             VANIS(2,IDIME,IELEX)=PARAM(1+NDIME+IDIME)             ! b0
            ENDDO
           ENDIF
          ENDIF
C
          IF(IFREN.EQ.56) THEN                       ! Myocardium model
           DO IDIME=1,3
            VANIS(1,IDIME,IELEX)=PARAM(1+IDIME)                    ! f0
            VANIS(2,IDIME,IELEX)=PARAM(4+IDIME)                    ! s0
           ENDDO
          ENDIF
C
          IF(IFREN.EQ.2.OR.IFREN.EQ.3) THEN          ! standard models
           DO IDIME=1,NDIME
            VANIS(1,IDIME,IELEX)=PARAM(1+IDIME)                    ! vx
            IF(NDIME.EQ.3) VANIS(2,IDIME,IELEX)=PARAM(4+IDIME)     ! vy
           ENDDO
           IF(NDIME.EQ.2) THEN
            VANIS(2,1,IELEX)=-VANIS(1,2,IELEX)
            VANIS(2,2,IELEX)= VANIS(1,1,IELEX)
           ENDIF
          ENDIF
C
          IF(NFIBN.GT.0) THEN
           IY=2
           IF(NDIME.EQ.3) IY=1
           DO I=1,NFIBN
            DO IDIME=1,NDIME
             IX=(I-1)*(NDIME+1)+IDIME
             VANIS(I,IDIME,IELEX)=PARAM(1+IX)                  ! N_theta
            ENDDO
            IX=I*(NDIME+1)
            VANIS(I,NDIME+IY,IELEX)=PARAM(1+IX)                ! f_theta
           ENDDO
          ENDIF
C
         ELSE
          CALL RUNEND('PROPANI: ISOTR. ELEMENT IN MATER._SYST._COORD.')
         ENDIF
        ENDDO
C
        CALL LISTEN('PROPANI',NPRIN,ITAPE)
        IF(WORDS(1).NE.'END_M')
     .   CALL RUNEND('PROPANI:END_MATERIAL_SYST._COORD. CARD NOT FOUND')
C
        WRITE(LURES,900)
        VATOL=1.0D-03
        DO IELEM=1,NELEM                   ! control
         LGRUP=MATNO(IELEM)                ! set number
         LMATS=INT(PROEL(1,LGRUP))         ! material number
         LTYPE=INT(PROEL(5,LGRUP))
         IF(LTYPE.EQ.30) THEN
          CALL IDENPR(PROPS(1,LMATS))      ! obtains material properties
          ISOTR=INT(PROPS(1,LMATS))        ! 0=isotropic; 1=orthotropic
          IF(ISOTR.EQ.1) THEN
           IF(IFREN.EQ.2.OR.IFREN.EQ.3.OR.
     .        IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57) THEN
            DO I=1,2                       ! unit vectors
             VANIM=VANIS(I,1,IELEM)*VANIS(I,1,IELEM)+
     .             VANIS(I,2,IELEM)*VANIS(I,2,IELEM)+
     .             VANIS(I,3,IELEM)*VANIS(I,3,IELEM)
             IF(DABS(VANIM-1.0D0).GT.VATOL) THEN
              WRITE(LURES,912) IELEM
              CALL RUNEND('PROPANI: ELEMENT WITH NO MATER_SYST._COORD.')
             ENDIF
             IF(IFREN.EQ.55.OR.IFREN.EQ.56.OR.IFREN.EQ.57) THEN
              WRITE(LURES,901) IELEM,I,(VANIS(I,IDIME,IELEM), IDIME=1,3)
             ENDIF
             IF(IFREN.EQ.2.OR.IFREN.EQ.3) THEN
              WRITE(LURES,902) IELEM,I,(VANIS(I,IDIME,IELEM),
     .                                                    IDIME=1,NDIME)
             ENDIF
            ENDDO
           ENDIF          ! ifren.eq.2.or.ifren.eq.3 ...
           IF(IFREN.EQ.2.OR.IFREN.EQ.3.OR.
     .        IFREN.EQ.56) THEN
            VANIM=VANIS(1,1,IELEM)*VANIS(2,1,IELEM)+   ! orthog. vectors
     .            VANIS(1,2,IELEM)*VANIS(2,2,IELEM)+
     .            VANIS(1,3,IELEM)*VANIS(2,3,IELEM)
            IF(DABS(VANIM).GT.VATOL) THEN
             WRITE(LURES,912) IELEM
             CALL RUNEND('PROPANI: f0 & s0 ARE NOT ORTHOGONAL VECTORS')
            ENDIF
           ENDIF          ! ifren.eq.2.or.ifren.eq.3.or.ifren.eq.56
           IF(NFIBN.GT.0) THEN
            VANIF=0.0D0
            DO I=1,NFIBN
             VANIM=VANIS(I,1,IELEM)*VANIS(I,1,IELEM)+
     .             VANIS(I,2,IELEM)*VANIS(I,2,IELEM)+
     .             VANIS(I,3,IELEM)*VANIS(I,3,IELEM)
             IF(DABS(VANIM-1.0D0).GT.VATOL) THEN
              WRITE(LURES,912) IELEM
              CALL RUNEND('PROPANI: ELEMENT WITH NO MATER_SYST._COORD.')
             ENDIF
             VANIF=VANIF+VANIS(I,4,IELEM)
            ENDDO
            IF(DABS(VANIF-1.0D0).GT.VATOL) THEN
             WRITE(LURES,912) IELEM
             CALL RUNEND('PROPANI: ELEMENT WITH WRONG VOLUM. FRACTIONS')
            ENDIF
           ENDIF          ! nfibn.gt.0
          ENDIF           ! isotr.eq.1
         ENDIF            ! ltype.eq.30
        ENDDO             ! ielem=1,nelem
       ELSE
        CALL RUNEND('INPPRO: NO MAT. SYSTEM OF COORD. FOR ANISOTROPY')
       ENDIF
      ELSE
       IF(NANIS.EQ.1)
     .  CALL RUNMEN('WARNING: NANIS=1 IS NOT NECESSARY (NO ANIS. MODEL')
      ENDIF               ! ianis.eq.1
C
  900 FORMAT(//,'MATERIAL SYSTEM OF COORDINATES FOR ANISOTROPY')
  901 FORMAT('ELEMENT=',I6,2X,'FIBER NUMBER=',I2,2X,
     .       'FIBER DIRECTION=',3E15.6)
  902 FORMAT('ELEMENT=',I6,2X,'AXIS NUMBER=',I2,2X,
     .       'AXIS DIRECTION=',3E15.6)
  912 FORMAT(//,'ERROR: ANIS. ELEMENT WITH NO MATERIAL_SYST._COORD.=',
     .       I6)
      END
