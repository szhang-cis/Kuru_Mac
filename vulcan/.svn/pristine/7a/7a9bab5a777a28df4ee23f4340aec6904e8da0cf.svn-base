      PROGRAM MLITMESH
C***********************************************************************
C
C**** THIS ROUTINE HAS BEEN BUILT FROM:
C     addworf.f, remfro.f & mlirt3d.f
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(MAXCP=5300, MDIMEF=3,              ! change if necessary
     .          MWORK=85000000, MMOVI=1, MPREC=1, MNCEL=3, MAXCE=5000 )
C
      DIMENSION ICONTF(35), X(MAXCP+6), Y(MAXCP+6),
     .          CCONT(MDIMEF+2+MMOVI*(MDIMEF+MPREC),MAXCP),
     .          IPONT(20),NPONT(21),FBC(20,3)
C
      CHARACTER*40 FILEIN1,FILEOUT1
      
      real(kind=8),allocatable :: WORK(:), LNODC(:,:)
      
      allocate(WORK(MWORK))
      allocate(LNODC(MNCEL+1,MAXCE))
C
C**** ROUTINE CONMOV
C
      KREMF=1
      NDIMEF=3             ! le MDIMEF !
      KMOVI=1              ! eq MMOVI !
      KPREC=1
      ISTEPF=1
      LURESF=83            ! fileout1
C
      write(6,*)'Input file 1 name: ?'
      read(5,*) filein1
      write(6,*)'Output file 1 name: ?'
      read(5,*) fileout1

      open(51,file=filein1)
      open(83,file=fileout1)

      write(6,*) 'Number of contour points & elements'
      read(5,*) ncopo,ncoel
      if(ncopo.gt.maxcp) then
       write(6,*) '**** ncopo gt maxcp ****',ncopo,maxcp
       stop
      endif
      if(ncoel.gt.maxce) then
       write(6,*) '**** ncoel gt maxce ****',ncoel,maxce
       stop
      endif
C
      MELCP=20
      NNCEL=3              ! le MNCEL !
      MARIS=3*MAXCE        ! rough estimation
      KPARF=1
C
c     IBUBL=0              ! 2D
c     NCOPN=NCOPO
c     IDGRE=2
C
      icontf(9)=0          !!
      IWORK=ICONTF( 9)-1
C
      IF(KREMF.EQ.1) THEN
c      IF(NDIMEF.EQ.2) THEN
c       NNCOP=MAXCP+6   ! NCOPO+2*IDELT, IDELT=0 OPEN, =3 CLOSED CONT.
c       NC=MAXCP        ! assumption
c       ID=IDGRE+1
c       ICONTF( 9)=ICONTF( 9)               ! x(NNCOP)
c       ICONTF(10)=ICONTF( 9)+NNCOP         ! y(NNCOP)
c       ICONTF(11)=ICONTF(10)+NNCOP         ! tn(NNCOP)
c       ICONTF(12)=ICONTF(11)+NNCOP         ! ratio(NNCOP,2)
c       ICONTF(13)=ICONTF(12)+2*NNCOP       ! te(NNCOP)
c       ICONTF(14)=ICONTF(13)+NNCOP         ! AN(3,nc)
c       ICONTF(15)=ICONTF(14)+3*nc          ! curv(NNCOP)
c       ICONTF(16)=ICONTF(15)+NNCOP         ! xnew(NNCOP)
c       ICONTF(17)=ICONTF(16)+NNCOP         ! ynew(NNCOP)
c       ICONTF(18)=ICONTF(17)+NNCOP         ! x(NNCOP)
c       ICONTF(19)=ICONTF(18)+NNCOP         ! y(NNCOP)
c       ICONTF(20)=ICONTF(19)+NNCOP         ! ts(NNCOP)
c       ICONTF(21)=ICONTF(20)+NNCOP         ! W(nc,nc)
c       ICONTF(22)=ICONTF(21)+nc*nc         ! P(idegre+1,nc)
c       ICONTF(23)=ICONTF(22)+id*nc         ! D(idegre+1,idegre+1)
c       ICONTF(24)=ICONTF(23)+id*id         ! aux1(idegre+1,nc)
c       ICONTF(25)=ICONTF(24)+id*nc         ! aux2(idegre+1,idegre+1)
c       ICONTF(26)=ICONTF(25)+id*id         ! aux3(idegre+1,idegre+1)
c       ICONTF(27)=ICONTF(26)+id*id         ! v(NNCOP+6)
c       ICONTF(28)=ICONTF(27)+NNCOP+6       ! s(NNCOP) num. of intervals
c       ICONTF(29)=ICONTF(28)+NNCOP         ! interpc(NNCOP)  curvature
c       ICONTF(30)=ICONTF(29)+(NNCOP+1)/2   ! interpm(NNCOP)  meshing
c       LCONTF    =ICONTF(30)+(NNCOP+1)/2
c      END IF
       IF(NDIMEF.EQ.3) THEN
        ICONTF(10)=ICONTF( 9)
        ICONTF(11)=ICONTF(10)                 ! xn(maxcp,ndimef)
        ICONTF(12)=ICONTF(11)+MAXCP*NDIMEF    ! nnube(maxcp,melcp*nncel)
        ICONTF(13)=ICONTF(12)+(MAXCP*MELCP*NNCEL+1)/2   ! ncpnu(maxcp,5)
        ICONTF(14)=ICONTF(13)+(MAXCP*5+1)/2   ! lelnu(maxcp,melcp)
        ICONTF(15)=ICONTF(14)+(MAXCP*MELCP+1)/2         ! krelm(maris,2)
        ICONTF(16)=ICONTF(15)+(MARIS*2+1)/2   ! krnoe(maris,2)
        ICONTF(17)=ICONTF(16)+(MARIS*2+1)/2   ! kg(maxcp,maxcp)
        ICONTF(18)=ICONTF(17)+(MAXCP*MAXCP+1)/2         ! kelar(maxce,3)
        ICONTF(19)=ICONTF(18)+(MAXCE*3+1)/2   ! xnewn(maris,ndimef)
        ICONTF(20)=ICONTF(19)+MARIS*NDIMEF    ! kindic(maris,2)
        ICONTF(21)=ICONTF(20)+(MARIS*2+1)/2   ! xpcp(melcp*2+2,3)
        ICONTF(22)=ICONTF(21)+(MELCP*2+2)*3   ! h(melcp*2+2)
        ICONTF(23)=ICONTF(22)+MELCP*2+2       ! psi(melcp*2+2)
        ICONTF(24)=ICONTF(23)+MELCP*2+2       ! eta(melcp*2+2)
        ICONTF(25)=ICONTF(24)+MELCP*2+2       ! p(melcp*2+2,6)
        ICONTF(26)=ICONTF(25)+(MELCP*2+2)*6   ! xnn(maxcp,2*ndimef)
        ICONTF(27)=ICONTF(26)+MAXCP*2*NDIMEF  ! lnodn(nncel+1,maxce)
        ICONTF(28)=ICONTF(27)+                ! lnodl(nncel+1,melcp)
     .                        ((NNCEL+1)*MAXCE+1)/2
        ICONTF(29)=ICONTF(28)+                ! lnodp(melcp+1)
     .                        ((NNCEL+1)*MELCP+1)/2
        ICONTF(30)=ICONTF(29)+(MELCP+2)/2     ! mcoel(melcp)
        LCONTF    =ICONTF(30)+(MELCP+1)/2
       END IF
      END IF
c     IF(LCONTF.GT.LWOR1F) LWOR1F=LCONTF
      if(lcontf.gt.mwork) then
       write(6,*) '**** lcontf gt mwork ****',lcontf,mwork
       stop
      endif
C
C**** READ CONTOUR MESH
C
C     CCONT(1:NDIMEF,IPOIN)           ! coordinates
C     CCONT(NDIMEF+2,IPOIN)           ! reference element to search
C     CCONT(NDIMEF+3,IPOIN)           ! prescribed x velocity
C     CCONT(NDIMEF+4,IPOIN)           ! prescribed y velocity
C     CCONT(NDIMEF+5,IPOIN)           ! prescribed z velocity
C     CCONT(NDIMEF+6,IPOIN)           ! prescribed pressure
C
      read(51,*)                      ! CONTOUR
      read(51,*)                      ! COORDINATES
      do i=1,ncopo
       read(51,*) ipoin,(CCONT(j,IPOIN), j=1,ndimef),
     .                  (CCONT(NDIMEF+1+k,IPOIN), k=1,5)
      end do
      read(51,*)                      ! END_COORDINATES
      read(51,*)                      ! ELEMENTS
      do i=1,ncoel
       read(51,*) ielem,(lnodc(j,ielem),j=1,nncel),lnodc(nncel+1,ielem)
      enddo
      read(51,*)                      ! END_ELEMENTS
      read(51,*)                      ! END_CONTOUR

  903 format(i8,3e15.6,i8,4i3)
C
C**** ASSUMED PARAMETERS (change if necessary!)
C
      armin=0.020d0    ! minimum edge length
      kcute=1          ! =0: cut all edges; =1: cut ii & b edges
      proga=1.3d0      ! maximum area relation between adjacent elements
C
      anglm=5.0d0      ! minimum angle between adjacent elements
      anglb=60.0d0     ! maximum angle between adjacent boundary edges
      weiml=0.05d0     ! weight factor for baricentric repositioning
C
      KREDU=2
C
C**** INTERFACE REMESHING
C
c     IF(NDIMEF.EQ.2)
c    . CALL MLIRT2D(CCONT,NCOPO,NCOEL,NDIMEF,KMOVI,KPREC,MAXCP,
c    .              IBUBL,NCOPN,IDGRE,X,Y,WORK,ICONTF,LURESF,
c    .              IWORK,KPARF)
      IF(NDIMEF.EQ.3)
     . CALL MLIRT3D(CCONT,NCOPO,NCOEL,NDIMEF,KMOVI,KPREC,MAXCP,
     .              IBUBL,      IDGRE,                LURESF,
     .                    KPARF,
     .              MAXCE,NNCEL,MELCP,MARIS,KREDU,
     .              LNODC,ISTEPF,ARMIN,KCUTE,PROGA,ANGLM,ANGLB,
     .              WEIML,
     .              WORK(ICONTF(11)-IWORK),WORK(ICONTF(12)-IWORK),
     .              WORK(ICONTF(13)-IWORK),WORK(ICONTF(14)-IWORK),
     .              WORK(ICONTF(15)-IWORK),WORK(ICONTF(16)-IWORK),
     .              WORK(ICONTF(17)-IWORK),WORK(ICONTF(18)-IWORK),
     .              WORK(ICONTF(19)-IWORK),WORK(ICONTF(20)-IWORK),
     .              WORK(ICONTF(21)-IWORK),WORK(ICONTF(22)-IWORK),
     .              WORK(ICONTF(23)-IWORK),WORK(ICONTF(24)-IWORK),
     .              WORK(ICONTF(25)-IWORK),WORK(ICONTF(26)-IWORK),
     .              WORK(ICONTF(27)-IWORK),WORK(ICONTF(28)-IWORK),
     .              WORK(ICONTF(29)-IWORK))
C
      STOP
      END
C***********************************************************************
      SUBROUTINE MLIRT3D(CCONT,NCOPO,NCOEL,NDIMEF,KMOVI,KPREC,MAXCP,
     .                   IBUBL,      IDGRE,                LURESF,
     .                         KPARF,
     .                   MAXCE,NNCEL,MELCP,MARIS,KREDU,
     .                   LNODC,ISTEPF,ARMIN,KCUTE,PROGA,ANGLM,ANGLB,
     .                   WEIML,
     .                   xn,nnube,ncpnu,lelnu,krelm,krnoe,
     .                   kg,kelar,xnewn,kindic,xpcp,h,
     .                   psi,eta,p,xnn,lnodn,lnodl,lnodp)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS THE FRONT REMESHING IN 3D
C
C***********************************************************************
c
c     Notes:
c
c     nnube dimension: ncopo x (nelem surrounding a contour point
c                      x number of nodes per element (3 for triangles))
C
C     LIST OF VARIABLES (WORK ARRAYS):
C     xn(maxcp,ndimef)                         (equals CCONT)
C     nnube(maxcp,melcp*nncel)
C     ncpnu(maxcp,5)
C     lelnu(maxcp,melcp)
C     krelm(maris,2)
C     krnoe(maris,2)
C     kg(maxcp,maxcp)
C     kelar(maxce,3)
C     xnewn(maris,ndimef)
C     kindic(maris,2)
C     xpcp(melcp*2+2,3)
C     h(melcp*2+2)
C     psi(melcp*2+2)
C     eta(melcp*2+2)
C     p(melcp*2+2,6)
C     xnn(maxcp,2*ndimef)
C     lnodn(nncel+1,maxce)
C     lnodl(nncel+1,melcp)
C     lnodp(melcp+1)
C
C     LIST OF PARAMETERS:
C     NCOPO: number of contour points (cp)
C     NCOEL: number of contour elements (ce)
C     MAXCP: maximum number of cp
C     MAXCE: maximum number of ce
C     IBUBL: not used
C     IDGRE: polinomial degree
C     KPARF: not used
C     NNCEL: number of nodes per contour element (=3)
C     MELCP: maximum number of surrounding elements to a given cp
C     MARIS: maximum number of edges
C     KREDU: =0 no cp reduction; =1 one cp reduction; =2 recursive red.
C
C
C     Nomenclature:
C
C     Interface internal point: surrounded by interface elements (ii)
C     Boundary internal point: surrounded by boundary elements (bi)
C     Boundary point: surrounded by interface & boundary elements (b)
C
C     Interface element: element of the (moving) interface (lnodc(4)=0)
C     Boundary element: element of the boundary (lnodc(4)=1)
C
C     Interface internal edge: shared by two interface elements
C     Boundary internal edge: shared by two boundary elements
C     Boundary edge: shared by an interface element & a boundary element
C
C
C     Steps after initial cloud identification:
C
C     1) Cut edges according to angle & area ratio criteria
C        Only for interface internal edges
C     2) Do not cut edges according to a minimum length
C        For all edges
C     3) Cut edges according to a maximum length
C        kcute=0: for all edges; kcute=1: for ii & b edges
C     4) Computes new coordinates at middle of cut edges, builds new
C        contour coordinate and connectivity arrays & redefines cloud
C     5) Improves quality by diagonal inversion according to angle,
C        quality, convexity and area ratio criteria
C        Only for interface and boundary internal edges
C     6) Check for reducing number of contour points according to two
C        criteria based on coplanarity and minimum edge length
C        (redefines cloud, improves quality by diagonal inversion &
C        apply baricentric repositioning)
C        Only for interface internal contour points
C     7) Apply patch-based baricentric repositioning
C        Only for interface internal contour points
C        Only for boundary internal contour points according to boundary
C        conditions (bc)
C     8) Apply edge-based baricentric repositioning according to bc
C        Only for boundary contour points
C
C***********************************************************************
      implicit double precision (a-h,o-z)
C
      DIMENSION CCONT(NDIMEF+2+KMOVI*(NDIMEF+KPREC),MAXCP)
      DIMENSION LNODC(NNCEL+1,MAXCE)
C
      dimension xn(maxcp,ndimef)           ! work array and input
      dimension nnube(maxcp,melcp*nncel),  ! work arrays
     .          ncpnu(maxcp,5),
     .          lelnu(maxcp,melcp),
     .          krelm(maris,2), krnoe(maris,2), kg(maxcp,maxcp),
     .          kelar(maxce,3), xnewn(maris,ndimef),
     .          kindic(maris,2), xpcp(melcp*2+2,3), h(melcp*2+2),
     .          psi(melcp*2+2), eta(melcp*2+2), p(melcp*2+2,6)
      dimension xnn(maxcp,2*ndimef),       ! work arrays and then output
     .          lnodn(nncel+1,maxce)
      dimension lnodl(nncel+1,melcp), lnodp(melcp+1)
C
      DIMENSION edge1(3), edge2(3), edge3(3), areai(3), areaj(3),
     .          t1(3), t2(3),
     .          paux(6,6), acoef(6), pauy(3,3),
     .          xpnew(3)
      dimension w(3),fv1(3)        ! auxiliars for eigenvalue RS rutines
      dimension icopl(4)
C
      TWOPI=6.283185307179586D0
C
C**** ASSUMED PARAMETERS
C
      armax=5.0d0*armin            ! maximun edge length
c
      cosam=dcos(anglm*twopi/360.0d0)
      cosab=dcos(anglb*twopi/360.0d0)
C
C**** ASSIGNS COORDINATES
C
      DO IDIME=1,NDIMEF
       DO ICOPO=1,NCOPO
        XN(ICOPO,IDIME)=CCONT(IDIME,ICOPO)
       END DO
      END DO
C
C**** GiD OUTPUT BEFORE REMESHING
C
c     CALL OUTGID(LNODC,XN,CCONT,
c    .            MAXCE,NNCEL,NDIMEF,KMOVI,KPREC,MAXCP,
c    .            NCOEN,NCOPN,ISTEPF,    1)
C
C**** CLOUD IDENTIFICATION:
C
C     - number elements including the contour point icopo
C     - contour points belonging to them
C
C     nnube includes the point number belonging to adjacent elements
C
C     nnube=cnn (Piotr's nomenclature)
C
      CALL CLOUDI(LNODC,NNUBE,NCPNU,LELNU,KRELM,KRNOE,   KG,KELAR,
     .            MAXCE,MAXCP,MELCP,MARIS,
     .            NCOEL,NCOPO,NNCEL,KARIS,LURESF)
C
C**** STEP 1: Identifies the internal edges to be cut according to angle
C             & area ratio criteria
C
      do iaris=1,karis
       kindic(iaris,1)=0               ! identification index
       kindic(iaris,2)=0               ! new node number (used later)
      end do
      do iaris=1,karis                 ! computes angle between normals
       ie1=-1                          ! to deal with non-conforming el
       ie2=-1
       if(krelm(iaris,2).ne.0) then
        ie1=lnodc(nncel+1,krelm(iaris,1))
        ie2=lnodc(nncel+1,krelm(iaris,2))
       end if
       if(ie1.eq.0.and.ie2.eq.0) then  ! only interface internal edges
        i1=lnodc(1,krelm(iaris,1))     ! first element sharing the edge
        i2=lnodc(2,krelm(iaris,1))
        i3=lnodc(3,krelm(iaris,1))
        j1=lnodc(1,krelm(iaris,2))     ! second element sharing the edge
        j2=lnodc(2,krelm(iaris,2))
        j3=lnodc(3,krelm(iaris,2))
        do j=1,3                       ! 3D
         edge1(j)=xn(i2,j)-xn(i1,j)
         edge2(j)=xn(i3,j)-xn(i1,j)
        end do
        call vecpro(3,edge1,edge2,areai)
        areaim=areai(1)*areai(1)+areai(2)*areai(2)+areai(3)*areai(3)
        do j=1,3                       ! 3D
         edge1(j)=xn(j2,j)-xn(j1,j)
         edge2(j)=xn(j3,j)-xn(j1,j)
        end do
        call vecpro(3,edge1,edge2,areaj)
        areajm=areaj(1)*areaj(1)+areaj(2)*areaj(2)+areaj(3)*areaj(3)
        do i=1,3                       ! normalization
         areai(i)=areai(i)/dsqrt(areaim)
         areaj(i)=areaj(i)/dsqrt(areajm)
        end do
        cosarea=areai(1)*areaj(1)+areai(2)*areaj(2)+areai(3)*areaj(3)
        if(cosarea.lt.cosam) then      ! non-(nearly) coplanar elements
         do j=1,3
          jaris=kelar(krelm(iaris,1),j)
          kindic(jaris,1)=1            ! edge cut
          jaris=kelar(krelm(iaris,2),j)
          kindic(jaris,1)=1            ! edge cut
         end do
        end if
C
        ratia=dsqrt(areaim/areajm)
        if(ratia.gt.proga.or.ratia.lt.(1.0d0/proga)) then
         if(ratia.gt.1.0d0) then                       ! area i > area j
          i1=krelm(iaris,1)
          kindic(kelar(i1,1),1)=1                      ! edge cut
          kindic(kelar(i1,2),1)=1                      ! edge cut
          kindic(kelar(i1,3),1)=1                      ! edge cut
         else                                          ! area j > area i
          i1=krelm(iaris,2)
          kindic(kelar(i1,1),1)=1                      ! edge cut
          kindic(kelar(i1,2),1)=1                      ! edge cut
          kindic(kelar(i1,3),1)=1                      ! edge cut
         end if
        end if
       end if            ! ie1=0 & ie2=0
      end do             ! karis
C
C**** Identifies the boundary edges to be cut according to their link
C     to internal edges (this is a Piotr's criterion to be revised!)
C
      do iaris=1,karis
       ie1=-1                          ! to deal with non-conforming el
       ie2=-1
       if(krelm(iaris,2).ne.0) then
        ie1=lnodc(nncel+1,krelm(iaris,1))
        ie2=lnodc(nncel+1,krelm(iaris,2))
       end if
c      if((ie1.eq.1.and.ie2.eq.0).or.(ie1.eq.0.and.ie2.eq.1)) then
c       i1=krelm(iaris,1)
c       if(kindic(kelar(i1,1),1).eq.1.or.kindic(kelar(i1,2),1).eq.1.or.
c    .     kindic(kelar(i1,3),1).eq.1) kindic(iaris,1)=1
c      end if
C
C**** STEPS 2 & 3: Identifies length of all edges
C
       do j=1,3                        ! 3D
        edge1(j)=xn(krnoe(iaris,2),j)-xn(krnoe(iaris,1),j)
       end do
       edgem=edge1(1)*edge1(1)+edge1(2)*edge1(2)+edge1(3)*edge1(3)
       edgem=dsqrt(edgem)
       if(edgem.lt.(2.d0*armin)) kindic(iaris,1)=0
       if(kcute.eq.0) then             ! cut all edges (ii, bi & b)
        if(edgem.gt.(2.d0*armax)) kindic(iaris,1)=1
       end if
       if(kcute.eq.1) then             ! cut ii & b edges
        if(ie1.eq.0.or.ie2.eq.0) then
         if(edgem.gt.(2.d0*armax)) kindic(iaris,1)=1
        end if
       end if
C
C**** STEP 4: Computes new coordinates at middle of cut edges
C
       do j=1,3
        xnewn(iaris,j)=0.0d0
       end do
       if(kindic(iaris,1).eq.1) then
        do j=1,3
         xnewn(iaris,j)=xn(krnoe(iaris,1),j)+
     .                 (xn(krnoe(iaris,2),j)-xn(krnoe(iaris,1),j))/2.0d0
        end do
C
C**** Computes new coordinates satisfying the least square distance
C     function and selects plane of projection at the added point
C     according to a normal vector computed as average of normal vectors
C     at elements that share the edge, t1 (suggested) and t2=n x t1
C

        go to 9999   ! to be revised (to take into account bound. cond.)

        if(ie1.eq.0.and.ie2.eq.0) then ! only interface internal edges
         ncpe=0
         do i=1,ncpnu(krnoe(iaris,1),1)
          i1=nnube(krnoe(iaris,1),i)
          if(ncpnu(i1,5).ne.2) ncpe=ncpe+1      ! ii or b contour points
         end do
c
         do i=1,ncpnu(krnoe(iaris,2),1)
          i2=nnube(krnoe(iaris,2),i)
          if(ncpnu(i2,5).ne.2) then
           kk=0
           do j=1,ncpnu(krnoe(iaris,1),1)
            i1=nnube(krnoe(iaris,1),j)
            if(ncpnu(i1,5).ne.2) then
             if(i2.eq.i1) kk=1         ! point of cloud 2 included in 1
            end if
           end do
           if(kk.eq.0) ncpe=ncpe+1
          end if
         end do
         if(ncpe.ge.6) then            ! quadratic least squares
          i1=lnodc(1,krelm(iaris,1))   ! first element sharing the edge
          i2=lnodc(2,krelm(iaris,1))
          i3=lnodc(3,krelm(iaris,1))
          j1=lnodc(1,krelm(iaris,2))   ! second element sharing the edge
          j2=lnodc(2,krelm(iaris,2))
          j3=lnodc(3,krelm(iaris,2))
          do j=1,3                     ! 3D
           edge1(j)=xn(i2,j)-xn(i1,j)
           edge2(j)=xn(i3,j)-xn(i1,j)
          end do
          call vecpro(3,edge1,edge2,areai)
          areaim=areai(1)*areai(1)+areai(2)*areai(2)+areai(3)*areai(3)
          do j=1,3                     ! 3D
           edge1(j)=xn(j2,j)-xn(j1,j)
           edge2(j)=xn(j3,j)-xn(j1,j)
          end do
          call vecpro(3,edge1,edge2,areaj)
          areajm=areaj(1)*areaj(1)+areaj(2)*areaj(2)+areaj(3)*areaj(3)
c
          iopta=2         ! options to compute tangent (better as input)
          if(iopta.eq.1) then          ! 1st option to compute tangent
           do i=1,3                    ! edge norm vector (average)
            areai(i)=areai(i)/dsqrt(areaim)
            areaj(i)=areaj(i)/dsqrt(areajm)
            areai(i)=(areai(i)+areaj(i))/2.d0
            if(i.eq.1) then
             amin=dabs(areai(i))
             i1=1
            elseif(i.gt.1.and.dabs(areai(i)).lt.amin) then
             amin=dabs(areai(i))
             i1=i
            end if
           end do
           if(i1.eq.1) then
            t1(1)=0.d0
            t1(2)=areai(3)
            t1(3)=-areai(2)
           elseif(i1.eq.2) then
            t1(1)=areai(3)
            t1(2)=0.d0
            t1(3)=-areai(1)
           elseif(i1.eq.3) then
            t1(1)=areai(2)
            t1(2)=-areai(1)
            t1(3)=0.d0
           end if
          end if
          if(iopta.eq.2) then          ! 2nd option to compute tangent
           do i=1,3                    ! edge norm vector (average)
            areai(i)=areai(i)/dsqrt(areaim)
            areaj(i)=areaj(i)/dsqrt(areajm)
            areai(i)=(areai(i)+areaj(i))/2.d0
            t1(i)=xn(krnoe(iaris,2),i)-xn(krnoe(iaris,1),i)       ! edge
           end do
          end if
c
          call vecpro(3,areai,t1,t2)
          tm=dsqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
          t1(1)=t1(1)/tm
          t1(2)=t1(2)/tm
          t1(3)=t1(3)/tm
          tm=dsqrt(t2(1)*t2(1)+t2(2)*t2(2)+t2(3)*t2(3))
          t2(1)=t2(1)/tm
          t2(2)=t2(2)/tm
          t2(3)=t2(3)/tm
c
          ie=0
          do i=1,ncpnu(krnoe(iaris,1),1)
           i1=nnube(krnoe(iaris,1),i)
           if(ncpnu(i1,5).ne.2) then            ! ii or b contour points
            ie=ie+1
            do j=1,3
             xpcp(ie,j)=xn(i1,j)-xnewn(iaris,j)  ! changes coord. origin
            end do
           end if
          end do
c
          do i=1,ncpnu(krnoe(iaris,2),1)
           i2=nnube(krnoe(iaris,2),i)
           if(ncpnu(i2,5).ne.2) then
            kk=0
            do j=1,ncpnu(krnoe(iaris,1),1)
             i1=nnube(krnoe(iaris,1),j)
             if(ncpnu(i1,5).ne.2) then
              if(i2.eq.i1) kk=1        ! point of cloud 2 included in 1
             end if
            end do
            if(kk.eq.0) then
             ie=ie+1                ! number of points of the cloud edge
             do j=1,3
              xpcp(ie,j)=xn(i2,j)-xnewn(iaris,j)
             end do
            end if
           end if
          end do
          if(ncpe.ne.ie) then
c          call runendf('error: ncpe ne ie (mlirt3d.f)')
           write(83,*) 'error: ncpe ne ie (mlirt3d.f)'       ! mlit-mesh
           stop
          end if
c
          do i=1,ie                 ! loop over points of the cloud edge
           h(i)=0.d0
           psi(i)=0.d0
           eta(i)=0.d0
           do j=1,3                 ! projections over n,t1,t2
            h(i)=h(i)+xpcp(i,j)*areai(j)
            psi(i)=psi(i)+xpcp(i,j)*t1(j)
            eta(i)=eta(i)+xpcp(i,j)*t2(j)
           end do
           p(i,1)=1.0d0             ! least squares coefficient matrix
           p(i,2)=psi(i)
           p(i,3)=eta(i)
           p(i,4)=psi(i)*psi(i)/2.0d0
           p(i,5)=psi(i)*eta(i)
           p(i,6)=eta(i)*eta(i)/2.0d0
          end do
          do i=1,6
           do j=1,6
            paux(i,j)=0.d0
            do k=1,ie
             paux(i,j)=paux(i,j)+p(k,i)*p(k,j)
            end do
           end do
          end do
          call invert(paux,6,6)
          do i=1,6
           acoef(i)=0.d0
           do j=1,6
            do k=1,ie
             acoef(i)=acoef(i)+paux(i,j)*p(k,j)*h(k)
            end do
           end do
          end do
          do j=1,3
           xnewn(iaris,j)=xnewn(iaris,j)+acoef(1)*areai(j)
          end do
         end if           ! ncpe ge 6
        end if            ! ie1=0 & ie2=0
C
        if((ie1.eq.1.and.ie2.eq.0).or.      ! only boundary edges
     .     (ie1.eq.0.and.ie2.eq.1)) then
         ik=1
         if(ie2.eq.1) ik=2
         i1=lnodc(1,krelm(iaris,ik))  ! boundary element
         i2=lnodc(2,krelm(iaris,ik))
         i3=lnodc(3,krelm(iaris,ik))
         do j=1,3                     ! 3D
          edge1(j)=xn(i2,j)-xn(i1,j)
          edge2(j)=xn(i3,j)-xn(i1,j)
         end do
         call vecpro(3,edge1,edge2,areai)
         areaim=areai(1)*areai(1)+areai(2)*areai(2)+areai(3)*areai(3)
         do i=1,3                     ! outward normal to boundary elem.
          t2(i)=areai(i)/dsqrt(areaim)
         end do
c
         i1=krnoe(iaris,1)
         i2=krnoe(iaris,2)
         do j=1,3                        ! 3D
          t1(j)=xn(i2,j)-xn(i1,j)
         end do
         tm=dsqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
         do j=1,3                        ! 3D
          t1(j)=t1(j)/tm
         end do
         call vecpro(3,t1,t2,areai)
c
         ib=2                            ! assumed non-smooth boundary
         do ie=1,ncpnu(krnoe(iaris,1),1)
          i2=nnube(krnoe(iaris,1),ie)
          if(i2.ne.krnoe(iaris,1).and.
     .       i2.ne.krnoe(iaris,2)) then  ! looks for point of other edge
           if(ncpnu(i2,5).eq.1) then     ! boundary point of cloud 1
            i1=krnoe(iaris,1)
            do j=1,3                     ! 3D
             t2(j)=xn(i2,j)-xn(i1,j)
            end do
            tm=dsqrt(t2(1)*t2(1)+t2(2)*t2(2)+t2(3)*t2(3))
            do j=1,3                     ! 3D
             t2(j)=t2(j)/tm
            end do
            cosarea=t1(1)*t2(1)+t1(2)*t2(2)+t1(3)*t2(3)
            if(cosarea.gt.cosab) then
             ib=ib+1                     ! smooth boundary
             do j=1,3
              xpcp(ib,j)=xn(i2,j)-xnewn(iaris,j)
             end do
            end if
           end if
          end if
         end do
         do i=1,ncpnu(krnoe(iaris,2),1)
          i2=nnube(krnoe(iaris,2),i) 
          if(i2.ne.krnoe(iaris,1).and.
     .       i2.ne.krnoe(iaris,2)) then  ! looks for point of other edge
           if(ncpnu(i2,5).eq.1) then     ! boundary point of cloud 2
            i1=krnoe(iaris,2)
            do j=1,3                     ! 3D
             t2(j)=xn(i2,j)-xn(i1,j)
            end do
            tm=dsqrt(t2(1)*t2(1)+t2(2)*t2(2)+t2(3)*t2(3))
            do j=1,3                     ! 3D
             t2(j)=t2(j)/tm
            end do
            cosarea=t1(1)*t2(1)+t1(2)*t2(2)+t1(3)*t2(3)
            if(cosarea.gt.cosab) then
             ib=ib+1                     ! smooth boundary
             do j=1,3
              xpcp(ib,j)=xn(i2,j)-xnewn(iaris,j)
             end do
            end if
           end if
          end if
         end do
         if(ib.gt.2) then                ! least squares point found
          do j=1,3
           xpcp(1,j)=xn(krnoe(iaris,1),j)-xnewn(iaris,j)
           xpcp(2,j)=xn(krnoe(iaris,2),j)-xnewn(iaris,j)
          end do
c
          do i=1,ib                 ! number of points of the cloud edge
           h(i)=0.d0
           psi(i)=0.d0
           do j=1,3                 ! projections over n,t1,t2
            h(i)=h(i)+xpcp(i,j)*areai(j)
            psi(i)=psi(i)+xpcp(i,j)*t1(j)
           end do
           p(i,1)=1.0d0             ! least squares coefficient matrix
           p(i,2)=psi(i)
           p(i,3)=psi(i)*psi(i)/2.0d0
          end do
          do i=1,3
           do j=1,3
            pauy(i,j)=0.d0
            do k=1,ib
             pauy(i,j)=pauy(i,j)+p(k,i)*p(k,j)
            end do
           end do
          end do
          call invert(pauy,3,3)
          do i=1,3
           acoef(i)=0.d0
           do j=1,3
            do k=1,ib
             acoef(i)=acoef(i)+pauy(i,j)*p(k,j)*h(k)
            end do
           end do
          end do
          do j=1,3
           xnewn(iaris,j)=xnewn(iaris,j)+acoef(1)*areai(j)
          end do
         end if          ! ib gt 2
        end if           ! ie1=1 & ie2=0 or ...

 9999   continue

       end if            ! kindic(iaris,1) eq 1
      end do             ! iaris=1,karis
C
C**** Builds new contour coordinates array
C
      do icopo=1,ncopo
       do j=1,3
        xnn(icopo,j)=xn(icopo,j)
        xnn(icopo,ndimef+j)=CCONT(NDIMEF+2+j,ICOPO)
       end do
      end do
      newn=0
      do iaris=1,karis
       if(kindic(iaris,1).ne.0) then
        newn=newn+1
        kindic(iaris,2)=ncopo+newn
        do j=1,3
         xnn(ncopo+newn,j)=xnewn(iaris,j)
         xnn(ncopo+newn,ndimef+j)=0.0d0                   ! assumed free
         if(int(ccont(ndimef+2+j,krnoe(iaris,1))).eq.1.and.
     .      int(ccont(ndimef+2+j,krnoe(iaris,2))).eq.1)
     .      xnn(ncopo+newn,ndimef+j)=1.0d0           ! assumed criterion
        end do
       end if
      end do
      ncopn=ncopo+newn
      if(ncopn.gt.maxcp) then
c      CALL RUNENDF('ERROR: NUMBER OF NEW CONTOUR POINTS GT MAXCP')
       write(83,*) 'ncopn gt maxcp - ncopn=',ncopn           ! mlit-mesh
       stop
      end if
C
C**** Builds new contour connectivity array
C
      newe=0
      do icoel=1,ncoel
       i1=kelar(icoel,1)
       i2=kelar(icoel,2)
       i3=kelar(icoel,3)
       j1=0
       if(kindic(i1,1).ne.0) j1=1
       j2=0
       if(kindic(i2,1).ne.0) j2=2
       j3=0
       if(kindic(i3,1).ne.0) j3=4
       jj=j1+j2+j3
       if(jj.eq.0) then
        newe=newe+1
        do inods=1,nncel
         lnodn(inods,newe)=lnodc(inods,icoel)
        end do
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
       end if
       if(jj.eq.1) then
        newe=newe+1
        lnodn(1,newe)=lnodc(1,icoel)
        lnodn(2,newe)=kindic(i1,2)
        lnodn(3,newe)=lnodc(3,icoel)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=kindic(i1,2)
        lnodn(2,newe)=lnodc(2,icoel)
        lnodn(3,newe)=lnodc(3,icoel)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
       end if
       if(jj.eq.2) then
        newe=newe+1
        lnodn(1,newe)=lnodc(1,icoel)
        lnodn(2,newe)=kindic(i2,2)
        lnodn(3,newe)=lnodc(3,icoel)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=lnodc(1,icoel)
        lnodn(2,newe)=lnodc(2,icoel)
        lnodn(3,newe)=kindic(i2,2)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
       end if
       if(jj.eq.3) then
        newe=newe+1
        lnodn(1,newe)=lnodc(1,icoel)
        lnodn(2,newe)=kindic(i1,2)
        lnodn(3,newe)=lnodc(3,icoel)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=kindic(i1,2)
        lnodn(2,newe)=lnodc(2,icoel)
        lnodn(3,newe)=kindic(i2,2)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=kindic(i1,2)
        lnodn(2,newe)=kindic(i2,2)
        lnodn(3,newe)=lnodc(3,icoel)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
       end if
       if(jj.eq.4) then
        newe=newe+1
        lnodn(1,newe)=lnodc(1,icoel)
        lnodn(2,newe)=lnodc(2,icoel)
        lnodn(3,newe)=kindic(i3,2)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=kindic(i3,2)
        lnodn(2,newe)=lnodc(2,icoel)
        lnodn(3,newe)=lnodc(3,icoel)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
       end if
       if(jj.eq.5) then
        newe=newe+1
        lnodn(1,newe)=lnodc(1,icoel)
        lnodn(2,newe)=kindic(i1,2)
        lnodn(3,newe)=kindic(i3,2)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=kindic(i1,2)
        lnodn(2,newe)=lnodc(2,icoel)
        lnodn(3,newe)=lnodc(3,icoel)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=kindic(i3,2)
        lnodn(2,newe)=kindic(i1,2)
        lnodn(3,newe)=lnodc(3,icoel)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
       end if
       if(jj.eq.6) then
        newe=newe+1
        lnodn(1,newe)=lnodc(1,icoel)
        lnodn(2,newe)=lnodc(2,icoel)
        lnodn(3,newe)=kindic(i3,2)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=kindic(i3,2)
        lnodn(2,newe)=lnodc(2,icoel)
        lnodn(3,newe)=kindic(i2,2)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=kindic(i3,2)
        lnodn(2,newe)=kindic(i2,2)
        lnodn(3,newe)=lnodc(3,icoel)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
       end if
       if(jj.eq.7) then
        newe=newe+1
        lnodn(1,newe)=lnodc(1,icoel)
        lnodn(2,newe)=kindic(i1,2)
        lnodn(3,newe)=kindic(i3,2)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=kindic(i1,2)
        lnodn(2,newe)=lnodc(2,icoel)
        lnodn(3,newe)=kindic(i2,2)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=kindic(i3,2)
        lnodn(2,newe)=kindic(i2,2)
        lnodn(3,newe)=lnodc(3,icoel)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
        newe=newe+1
        lnodn(1,newe)=kindic(i1,2)
        lnodn(2,newe)=kindic(i2,2)
        lnodn(3,newe)=kindic(i3,2)
        lnodn(nncel+1,newe)=lnodc(nncel+1,icoel)
       end if
      end do                ! icoel=1,ncoel
      ncoen=newe
      if(ncoen.gt.maxce) then
c      CALL RUNENDF('ERROR: NUMBER OF NEW CONTOUR ELEM. GT MAXCE')
       write(83,*) 'ncoen gt maxce - ncoen=',ncoen           ! mlit-mesh
       stop
      end if
C
C**** REDEFINES CLOUD
C
      CALL CLOUDI(LNODN,NNUBE,NCPNU,LELNU,KRELM,KRNOE,   KG,KELAR,
     .            MAXCE,MAXCP,MELCP,MARIS,
     .            NCOEN,NCOPN,NNCEL,KARIS,LURESF)
C
C**** STEP 5: IMPROVES QUALITY BY DIAGONAL INVERSION
C
      CALL DIAGON(LNODN,NNUBE,NCPNU,LELNU,KRELM,KRNOE,KELAR,
     .            XNN,COSAM,
     .            MAXCE,MAXCP,MELCP,MARIS,
     .            NDIMEF,NNCEL,KARIS,LURESF)
C
C**** STEP 6: CHECK FOR REDUCING NUMBER OF CONTOUR POINTS
C
C     Notes:
C     kredu=1: collapse elements according to coplanarity
C              (to be revised!)
C     kredu=2: collapse edge nodes according to minimum length edge and
C              element coplanarity
C
      if(kredu.eq.1) then
       ncop1=ncopn
       ncop2=0                         ! always reduce for first time
       kredx=1
       do while ((ncop1-ncop2).gt.0.and.kredx.eq.1) 
C
        do icopo=1,ncopn
         nelcp=ncpnu(icopo,2)
         if(ncpnu(icopo,5).eq.0) then  ! assumed criterion (ii cp)
          ku=lelnu(icopo,1)            ! first element of the cloud
          i1=lnodn(1,ku)               ! normal to this element
          i2=lnodn(2,ku)
          i3=lnodn(3,ku)
          do j=1,3                     ! 3D
           edge1(j)=xnn(i2,j)-xnn(i1,j)
           edge2(j)=xnn(i3,j)-xnn(i1,j)
          end do
          call vecpro(3,edge1,edge2,areai)
          areaim=areai(1)*areai(1)+areai(2)*areai(2)+areai(3)*areai(3)
          areami=areaim                ! minimum area
          do i=1,3                     ! normalization
           areai(i)=areai(i)/dsqrt(areaim)
          end do
          kpr=1
          if(nelcp.gt.1) then
           do ielcp=2,nelcp            ! normal to the other elements
            ku=lelnu(icopo,ielcp)
            i1=lnodn(1,ku)
            i2=lnodn(2,ku)
            i3=lnodn(3,ku)
            do j=1,3                   ! 3D
             edge1(j)=xnn(i2,j)-xnn(i1,j)
             edge2(j)=xnn(i3,j)-xnn(i1,j)
            end do
            call vecpro(3,edge1,edge2,areaj)
            areajm=areaj(1)*areaj(1)+areaj(2)*areaj(2)+areaj(3)*areaj(3)
            if(areajm.lt.areami) areami=areajm
            do i=1,3                   ! normalization
             areaj(i)=areaj(i)/dsqrt(areajm)
            end do
            cosarea=areai(1)*areaj(1)+areai(2)*areaj(2)+
     .              areai(3)*areaj(3)
            if(cosarea.gt.cosam) kpr=kpr+1    ! nearly coplanar elements
           end do   ! ielcp
          end if
          if(kpr.eq.nelcp) ncpnu(icopo,3)=1          ! flag to remove cp
         end if
        end do     ! icopo
c
        do icopo=1,ncopn ! only remove points not belonging to the cloud
         if(ncpnu(icopo,3).eq.1) then
          ncp=ncpnu(icopo,1)
          do i=1,ncp
           j=nnube(icopo,i)
           if(j.ne.icopo.and.ncpnu(j,3).eq.1) ncpnu(j,3)=0 ! not removed
          end do
         end if
        end do     ! icopo
c
        do icopo=1,ncopn                               ! redefines lnodn
         if(ncpnu(icopo,3).eq.1) then
          ncp=ncpnu(icopo,1)
          nelcp=ncpnu(icopo,2)
          if((ncp-3).ge.2.or.((ncp-3).eq.1.and.(ncp-nelcp).eq.1)) then
           do ielcp=1,nelcp
            ku=lelnu(icopo,ielcp)
            do incel=1,nncel
             if(lnodn(incel,ku).eq.icopo) ix=incel
            end do
            lnodl(1,ielcp)=lnodn(ix,ku) ! always icopo in first position
            if(ix.eq.1) then
             lnodl(2,ielcp)=lnodn(2,ku)
             lnodl(3,ielcp)=lnodn(3,ku)
            end if
            if(ix.eq.2) then
             lnodl(2,ielcp)=lnodn(3,ku)
             lnodl(3,ielcp)=lnodn(1,ku)
            end if
            if(ix.eq.3) then
             lnodl(2,ielcp)=lnodn(1,ku)
             lnodl(3,ielcp)=lnodn(2,ku)
            end if
            lnodl(nncel+1,ielcp)=    ! =0 according to assumed criterion
     .      lnodn(nncel+1,ku)
           end do
c
           do ielcp=1,nelcp          ! builds polygon surrounding icopo
            jcp=2                    ! valid for boundary and internal
            lnodp(1)=lnodl(2,ielcp)  ! points
            lnodp(2)=lnodl(3,ielcp)
            do icp=1,ncp-3
             do jelcp=1,nelcp
              if(ielcp.ne.jelcp) then
               if(lnodl(2,jelcp).eq.lnodp(icp+1)) then
                lnodp(icp+2)=lnodl(3,jelcp)
                jcp=jcp+1            ! number of points of the polygon
               end if
              end if
             end do
            end do
            if(jcp.eq.(ncp-1)) go to 10
           end do
   10      continue
c
           if(jcp.ne.(ncp-1)) then                    ! obvious controls
c           call runendf('error jcp ne ncp-1')
            write(83,*) 'jcp.ne.(ncp-1)'                     ! mlit-mesh
            stop
           end if
           if((ncp-nelcp).eq.1.and.ielcp.ne.1) then
c           call runendf('error en punto interior')
            write(83,*) '(ncp-nelcp).eq.1.and.ielcp.ne.1'    ! mlit-mesh
            stop
           end if
c
c mejorar la programacion de esta parte! (e.g., poner un do icp=1,ncp-1)
c
           kpr=1
   12      do ielcp=1,ncp-3          ! checks convexity of the new patch
            lnodl(1,ielcp)=lnodp(1)
            do incel=2,nncel
             lnodl(incel,ielcp)=lnodp(incel+ielcp-1)
            end do
           end do
           i1=lnodl(1,1)             ! normal to this element
           i2=lnodl(2,1)
           i3=lnodl(3,1)
           do j=1,3                  ! 3D
            edge1(j)=xnn(i2,j)-xnn(i1,j)
            edge2(j)=xnn(i3,j)-xnn(i1,j)
           end do
           call vecpro(3,edge1,edge2,areai)
           areaim=areai(1)*areai(1)+areai(2)*areai(2)+areai(3)*areai(3)
           if(areaim.lt.1.0d-03*areami)
     .      go to 13                 ! does not remove aligned nodes
           do i=1,3                  ! normalization
            areai(i)=areai(i)/dsqrt(areaim)
           end do
           do ielcp=2,ncp-3          ! normal to the other elements
            i1=lnodl(1,ielcp)
            i2=lnodl(2,ielcp)
            i3=lnodl(3,ielcp)
            do j=1,3                 ! 3D
             edge1(j)=xnn(i2,j)-xnn(i1,j)
             edge2(j)=xnn(i3,j)-xnn(i1,j)
            end do
            call vecpro(3,edge1,edge2,areaj)
            areajm=areaj(1)*areaj(1)+areaj(2)*areaj(2)+areaj(3)*areaj(3)
            if(areajm.lt.1.0d-03*areami)
     .       go to 13                ! does not remove aligned nodes
            do i=1,3                 ! normalization
             areaj(i)=areaj(i)/dsqrt(areajm)
            end do
            proes=areai(1)*areaj(1)+areai(2)*areaj(2)+areai(3)*areaj(3)
            if(proes.lt.0.0d0)
     .       go to 13                ! avoids patch concavity
           end do    ! ielcp
           go to 14                  ! this node is going to be reduced
   13      kpr=kpr+1
           if(kpr.le.(ncp-1)) then   ! try another polygon
            iaux=lnodp(1)            ! change polygon node ordering
            do icp=2,ncp-1
             lnodp(icp-1)=lnodp(icp)
            end do
            lnodp(ncp-1)=iaux
            go to 12
           else
            ncpnu(icopo,3)=0         ! this node is not reduced
            go to 11
           end if
   14      continue
c
           do ielcp=1,nelcp
            ku=lelnu(icopo,ielcp)
            if(ielcp.le.(ncp-3)) then
             lnodn(1,ku)=lnodp(1)
             do incel=2,nncel
              lnodn(incel,ku)=lnodp(incel+ielcp-1)
             end do
             lnodn(nncel+1,ku)=      ! =0 according to assumed criterion
     .       lnodl(nncel+1,ielcp)
            else
             lnodn(1,ku)=0           ! flag to remove the element
            end if
           end do
          end if              ! ncp-3 ge 2 or ...
         end if               ! ncpnu(*,3)=1
   11    continue
        end do                ! icopo
c
        ncopx=ncopn
        do icopo=ncopn,1,-1
         if(ncpnu(icopo,3).eq.1) then
          do jcopo=icopo,ncopx
           do j=1,2*ndimef                 ! coordinates + prescriptions
            xnn(jcopo,j)=xnn(jcopo+1,j)
           end do
          end do
          ncopx=ncopx-1
         end if
        end do
        do icopo=1,ncopn
         ncpnu(icopo,4)=icopo
        end do
        do icopo=1,ncopn
         if(ncpnu(icopo,3).eq.1) then
          do jcopo=icopo,ncopn
           ncpnu(jcopo,4)=ncpnu(jcopo,4)-1
          end do
         end if
        end do
c
        ncoex=ncoen
        do icoel=ncoen,1,-1
         if(lnodn(1,icoel).eq.0) then    ! remove element
          do jcoel=icoel,ncoex
           do j=1,nncel
            lnodn(j,jcoel)=lnodn(j,jcoel+1)
           end do
           lnodn(nncel+1,jcoel)=lnodn(nncel+1,jcoel+1)
          end do
          ncoex=ncoex-1
         end if
        end do
        do icoel=1,ncoex
         do incel=1,nncel
          j=lnodn(incel,icoel)
          if(ncpnu(j,3).eq.0) lnodn(incel,icoel)=ncpnu(j,4)
         end do
        end do
c
        write(luresf,*) 'ncopn,ncopx=',ncopn,ncopx
        write(luresf,*) 'ncoen,ncoex=',ncoen,ncoex
c
        ncop1=ncopn                     ! useful for recursive reduction
        ncop2=ncopx
c
        ncopn=ncopx
        ncoen=ncoex
C
C**** REDEFINES CLOUD
C
        CALL CLOUDI(LNODN,NNUBE,NCPNU,LELNU,KRELM,KRNOE,   KG,KELAR,
     .              MAXCE,MAXCP,MELCP,MARIS,
     .              NCOEN,NCOPN,NNCEL,KARIS,LURESF)
C
C**** IMPROVES QUALITY BY DIAGONAL INVERSION
C
        CALL DIAGON(LNODN,NNUBE,NCPNU,LELNU,KRELM,KRNOE,KELAR,
     .              XNN,COSAM,
     .              MAXCE,MAXCP,MELCP,MARIS,
     .              NDIMEF,NNCEL,KARIS,LURESF)
C
C**** PERFORMS PATCH-BASED & EDGE-BASED BARICENTRIC REPOSITIONING
C
        CALL REPOSB(LNODN,NNUBE,NCPNU,LELNU,KRELM,KRNOE,
     .              XNN,COSAM,WEIML,
     .              NCOPN,
     .              MAXCE,MAXCP,MELCP,MARIS,
     .              NDIMEF,NNCEL,KARIS,LURESF)
C
        if(kredu.eq.1) kredx=0
       end do     ! while
      end if      ! kredu=1
C
      if(kredu.eq.2) then
       do iaris=1,karis
        ie1=-1                         ! to deal with non-conforming el
        ie2=-1
        if(krelm(iaris,2).ne.0) then
         ie1=lnodn(nncel+1,krelm(iaris,1))
         ie2=lnodn(nncel+1,krelm(iaris,2))
        end if
        ixa=0
        if((ie1.eq.0.and.ie2.eq.0).or.                   ! ii edge
     .     (ie1.eq.1.and.ie2.eq.1.and.kcute.eq.0)) then  ! bi edge
         if((ncpnu(krnoe(iaris,1),5).eq.0.and.           ! both ii nodes
     .       ncpnu(krnoe(iaris,2),5).eq.0).or.
     .      (ncpnu(krnoe(iaris,1),5).eq.2.and.           ! both bi nodes
     .       ncpnu(krnoe(iaris,2),5).eq.2)) ixa=1
        end if
        if((ie1.eq.1.and.ie2.eq.0).or.                   ! b edge
     .      (ie1.eq.0.and.ie2.eq.1)) then
         ixaa=0
         do j=1,ndimef
          if(int(xnn(krnoe(iaris,1),ndimef+j)).eq.
     .       int(xnn(krnoe(iaris,2),ndimef+j))) ixaa=ixaa+1
         end do
         if(ixaa.eq.ndimef) ixa=2        ! assumption: same bc of both n
        end if
        if(ixa.gt.0) then
         do j=1,3                        ! 3D
          edge1(j)=xnn(krnoe(iaris,2),j)-xnn(krnoe(iaris,1),j)
         end do
         edgem=edge1(1)*edge1(1)+edge1(2)*edge1(2)+edge1(3)*edge1(3)
         edgem=dsqrt(edgem)
c
         do ix=1,2                       ! evaluates coplanarity
          icopl(ix)=0                    ! computes normals
          nelcp=ncpnu(krnoe(iaris,ix),2) ! elements of the cloud
          nelcpx=0
          kpr=1
          do ielcp=1,nelcp               ! normal of the elements
           ku=lelnu(krnoe(iaris,ix),ielcp)
           ixb=0
           if((ixa.eq.1).or.
     .        (ixa.eq.2.and.lnodn(nncel+1,ku).eq.0)) then
            nelcpx=nelcpx+1              ! assumption: only ii elements
            ixb=1
           end if
           if(ixb.eq.1) then             ! normal to this element
            i1=lnodn(1,ku)
            i2=lnodn(2,ku)
            i3=lnodn(3,ku)
            do j=1,3                     ! 3D
             edge1(j)=xnn(i2,j)-xnn(i1,j)
             edge2(j)=xnn(i3,j)-xnn(i1,j)
            end do
            if(ielcp.eq.1) then
             call vecpro(3,edge1,edge2,areai)
             areaim=areai(1)*areai(1)+areai(2)*areai(2)+
     .                                                 areai(3)*areai(3)
             do i=1,3                    ! normalization
              areai(i)=areai(i)/dsqrt(areaim)
             end do
            else
             call vecpro(3,edge1,edge2,areaj)
             areajm=areaj(1)*areaj(1)+areaj(2)*areaj(2)+
     .                                                 areaj(3)*areaj(3)
             do i=1,3                    ! normalization
              areaj(i)=areaj(i)/dsqrt(areajm)
             end do
             cosarea=areai(1)*areaj(1)+areai(2)*areaj(2)+
     .                                                 areai(3)*areaj(3)
             if(cosarea.gt.cosam) kpr=kpr+1  ! nearly coplanar elements
            end if
           end if
          end do   ! ielcp
          if(kpr.eq.nelcpx) icopl(ix)=1
         end do
c
         do ix=1,2                       ! evaluates copl. of collap. n
          do j=1,ndimef                  ! stores coordinates
           xpnew(j)=xnn(krnoe(iaris,ix),j)
           xnn(krnoe(iaris,ix),j)=(xnn(krnoe(iaris,1),j)+
     .                             xnn(krnoe(iaris,2),j))/2.0d0
          end do
          icopl(ix+2)=0                  ! computes normals
          nelcp=ncpnu(krnoe(iaris,ix),2) ! elements of the cloud
          nelcpx=0
          kpr=1
          do ielcp=1,nelcp               ! normal of the elements
           ku=lelnu(krnoe(iaris,ix),ielcp)
           ixb=0
           if((ixa.eq.1).or.
     .        (ixa.eq.2.and.lnodn(nncel+1,ku).eq.0)) then
            nelcpx=nelcpx+1              ! assumption: only ii elements
            ixb=1
           end if
           if(ixb.eq.1) then             ! normal to this element
            i1=lnodn(1,ku)
            i2=lnodn(2,ku)
            i3=lnodn(3,ku)
            do j=1,3                     ! 3D
             edge1(j)=xnn(i2,j)-xnn(i1,j)
             edge2(j)=xnn(i3,j)-xnn(i1,j)
            end do
            if(ielcp.eq.1) then
             call vecpro(3,edge1,edge2,areai)
             areaim=areai(1)*areai(1)+areai(2)*areai(2)+
     .                                                 areai(3)*areai(3)
             do i=1,3                    ! normalization
              areai(i)=areai(i)/dsqrt(areaim)
             end do
            else
             call vecpro(3,edge1,edge2,areaj)
             areajm=areaj(1)*areaj(1)+areaj(2)*areaj(2)+
     .                                                 areaj(3)*areaj(3)
             do i=1,3                    ! normalization
              areaj(i)=areaj(i)/dsqrt(areajm)
             end do
             cosarea=areai(1)*areaj(1)+areai(2)*areaj(2)+
     .                                                 areai(3)*areaj(3)
             if(cosarea.gt.cosam) kpr=kpr+1  ! nearly coplanar elements
            end if
           end if
          end do   ! ielcp
          if(kpr.eq.nelcpx) icopl(ix+2)=1
          do j=1,ndimef                  ! restores coordinates
           xnn(krnoe(iaris,ix),j)=xpnew(j)
          end do
         end do
c
         if((edgem.lt.0.75d0*armin).or.                  ! 2 assumptions
     .      (edgem.lt.armax.and.
     .       icopl(1).eq.1.and.icopl(2).eq.1.and.
     .       icopl(3).eq.1.and.icopl(4).eq.1)) then      ! collapse node
          ix1=krnoe(iaris,1)          ! searches for lower numbering
          ix2=krnoe(iaris,2)
          if(ix1.gt.ix2) then
           ix1=krnoe(iaris,2)         ! lower numbering
           ix2=krnoe(iaris,1)         ! greater numbering
          end if
          do icoel=1,ncoen
           do incel=1,nncel
            if(lnodn(incel,icoel).eq.ix2) lnodn(incel,icoel)=ix1
            if(lnodn(incel,icoel).gt.ix2) lnodn(incel,icoel)=
     .                                              lnodn(incel,icoel)-1
           end do
          end do
c
          ncoex=ncoen                 ! reorders connectivity
          do icoel=ncoen,1,-1
           if(krelm(iaris,1).eq.icoel.or.krelm(iaris,2).eq.icoel) then
            if(icoel.ne.ncoex) then
             do jcoel=icoel+1,ncoex
              do incel=1,nncel+1
               lnodn(incel,jcoel-1)=lnodn(incel,jcoel)
              end do
             end do
            end if
            ncoex=ncoex-1
           end if
          end do
          if(ncoex.ne.(ncoen-2)) then
c          call runendf('error: ncoex ne ncoen-2')
           write(83,*) 'error: ncoex ne ncoen-2'            ! mlit-mesh
           stop
          end if
          ncoen=ncoex
c
          do j=1,ndimef               ! non-collapsed node coordinates
           xnn(ix1,j)=(xnn(ix1,j)+xnn(ix2,j))/2.0d0
          end do
          ncopx=ncopn                 ! reorders node numbering
          do icopn=ncopn,1,-1
           if(ix2.eq.icopn) then
            if(icopn.ne.ncopx) then
             do jcopn=icopn+1,ncopx
              do j=1,2*ndimef         ! coordinates + prescriptions
               xnn(jcopn-1,j)=xnn(jcopn,j)
              end do
             end do
            end if
            ncopx=ncopx-1
           end if
          end do
          if(ncopx.ne.(ncopn-1)) then
c          call runendf('error: ncopx ne ncopn-1')
           write(83,*) 'error: ncopx ne ncopn-1'             ! mlit-mesh
           stop
          end if
          ncopn=ncopx
C
C**** REDEFINES CLOUD
C
          CALL CLOUDI(LNODN,NNUBE,NCPNU,LELNU,KRELM,KRNOE,   KG,KELAR,
     .                MAXCE,MAXCP,MELCP,MARIS,
     .                NCOEN,NCOPN,NNCEL,KARIS,LURESF)
C
C**** IMPROVES QUALITY BY DIAGONAL INVERSION
C
          CALL DIAGON(LNODN,NNUBE,NCPNU,LELNU,KRELM,KRNOE,KELAR,
     .                XNN,COSAM,
     .                MAXCE,MAXCP,MELCP,MARIS,
     .                NDIMEF,NNCEL,KARIS,LURESF)
C
C**** PERFORMS PATCH-BASED & EDGE-BASED BARICENTRIC REPOSITIONING
C
          CALL REPOSB(LNODN,NNUBE,NCPNU,LELNU,KRELM,KRNOE,
     .                XNN,COSAM,WEIML,
     .                NCOPN,
     .                MAXCE,MAXCP,MELCP,MARIS,
     .                NDIMEF,NNCEL,KARIS,LURESF)
         end if   ! 2 assumptions
        end if    ! ixa gt 0
       end do     ! iaris=1,karis
      end if      ! kredu=2
C
C**** STEPS 7 & 8: PATCH-BASED & EDGE-BASED BARICENTRIC REPOSITIONING
C
      CALL REPOSB(LNODN,NNUBE,NCPNU,LELNU,KRELM,KRNOE,
     .            XNN,COSAM,WEIML,
     .            NCOPN,
     .            MAXCE,MAXCP,MELCP,MARIS,
     .            NDIMEF,NNCEL,KARIS,LURESF)
C
C**** ASSIGNS NEW NUMBER OF CP & CE, COORDINATES AND CONNECTIVITIES
C
      NCOPO=NCOPN
      NCOEL=NCOEN
      DO ICOPO=1,NCOPO
       DO IDIME=1,NDIMEF                ! coordinates + prescriptions
        CCONT(IDIME,ICOPO)=XNN(ICOPO,IDIME)
        CCONT(IDIME+NDIMEF+2,ICOPO)=XNN(ICOPO,IDIME+NDIMEF)
       END DO
       CCONT(NDIMEF+1,ICOPO)=0.0D0      ! inside the domain index
       CCONT(NDIMEF+2,ICOPO)=1.0D0      ! reference element to search
       CCONT(NDIMEF+6,ICOPO)=0.0D0      ! fix p      to be revised!
      END DO
      DO INCEL=1,NNCEL+1                ! connectivity + contour index
       DO ICOEL=1,NCOEL
        LNODC(INCEL,ICOEL)=LNODN(INCEL,ICOEL)
       END DO
      END DO
C
C**** GiD OUTPUT AFTER REMESHING
C
      CALL OUTGID(LNODN,XNN,CCONT,
     .            MAXCE,NNCEL,NDIMEF,KMOVI,KPREC,MAXCP,
     .            NCOEN,NCOPN,ISTEPF,    1)
C
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE CLOUDI(LNODC,NNUBE,NCPNU,LELNU,KRELM,KRNOE,   KG,KELAR,
     .                  MAXCE,MAXCP,MELCP,MARIS,
     .                  NCOEL,NCOPO,NNCEL,KARIS,LURESF)
C***********************************************************************
C
C**** CLOUD IDENTIFICATION:
C
C     - number elements including the contour point icopo
C     - contour points belonging to them
C
C     nnube includes the point number belonging to adjacent elements
C
C     nnube=cnn (Piotr's nomenclature)
C
C***********************************************************************
      implicit double precision (a-h,o-z)
C
      dimension lnodc(nncel+1,maxce)
      dimension nnube(maxcp,melcp*nncel),                  ! work arrays
     .          ncpnu(maxcp,5),
     .          lelnu(maxcp,melcp),
     .          krelm(maris,2), krnoe(maris,2), kg(maxcp,maxcp),
     .          kelar(maxce,3)
C
      do icopo=1,ncopo
       nelcp=0               ! number of elements adjacent to each icopo
       jcopn=0               ! number of nodes adjacent to each icopo
       do icoel=1,ncoel
        kpert=0
        do inods=1,nncel
         if(lnodc(inods,icoel).eq.icopo) kpert=1
        end do
        if(kpert.eq.1) then
         nelcp=nelcp+1
         if(nelcp.gt.melcp) then
          write(luresf,*)'numb.of elements adjacent to',icopo,'.gt.max'
          stop
         end if
         lelnu(icopo,nelcp)=icoel
         ncpnu(icopo,2)=nelcp
         do inods=1,nncel
          jcopn=jcopn+1
          nnube(icopo,jcopn)=lnodc(inods,icoel)
         end do
        end if
       end do
c
c**** Points number are duplicated in nnube, this duplication needs
c     to be eliminated
c     nnube = ncopo x number of points belonging to each icopo cloud
c
 20    do icopn=1,jcopn
        do kcopn=icopn+1,jcopn
         if(nnube(icopo,icopn).eq.nnube(icopo,kcopn)) then
          do lcopn=kcopn+1,jcopn
           nnube(icopo,lcopn-1)=nnube(icopo,lcopn)
          end do
          jcopn=jcopn-1
          go to 20
         end if
        end do
       end do
c
c**** Nomenclature:
c
c     ncpnu(1), number of points belonging to each icopo cloud:
c     - equal to the numer of triangles adjacent to icopo + 1, if icopo
c       does not belong to a boundary of the domain.
c     - equal to the numer of triangles adjacent to icopo + 2, if icopo
c       belongs to a boundary of the domain.
c
c     ncpnu(2), number of elements belonging to each icopo cloud
c
c     ncpnu(3), flag to remove contour point (1=to be removed)
c
c     ncpnu(4), new contour point number
c
c     ncpnu(5), index for (i & b) internal or boundary contour point
c               =0 interface internal (ii)
c               =1 boundary (b)
c               =2 boundary internal (bi)
c
c     ncpnu=ncnn (Piotr)
c
       ncpnu(icopo,1)=jcopn
c
      end do        ! icopo
C
      do icopo=1,ncopo
       ncpnu(icopo,3)=0       ! initialization
       ncpnu(icopo,4)=0
       ncpnu(icopo,5)=0
c
       nelcp=ncpnu(icopo,2)
       ix0=0
       ix1=0
       do ielcp=1,nelcp
        ku=lelnu(icopo,ielcp)
        ieu=lnodc(nncel+1,ku)
        if(ieu.eq.0) ix0=ix0+1
        if(ieu.eq.1) ix1=ix1+1
       end do
       if(ix0.ne.nelcp.and.ix1.ne.nelcp) ncpnu(icopo,5)=1
       if(ix1.eq.nelcp) ncpnu(icopo,5)=2
      end do
c
c**** Computes:
c
c     kg: auxiliar matrix
c     krelm: element numbers (2) that share each edge
c     krnoe: node numbers (2) that define each edge
c
c     do icopo=1,ncopo                 ! old form (to be removed!)
c      do jcopo=1,ncopo
c       kg(icopo,jcopo)=0
c      end do
c     end do
c     do icoel=1,ncoel
c      kg(lnodc(1,icoel),lnodc(2,icoel))=icoel
c      kg(lnodc(2,icoel),lnodc(3,icoel))=icoel
c      kg(lnodc(3,icoel),lnodc(1,icoel))=icoel
c     end do
C
c     karis=0
c     do icopo=1,ncopo-1
c      do jcopo=icopo+1,ncopo
c       if(kg(icopo,jcopo).ne.0.or.kg(jcopo,icopo).ne.0) then
c        karis=karis+1
c        krnoe(karis,1)=icopo
c        krnoe(karis,2)=jcopo
c        ii=0
c        if(kg(icopo,jcopo).ne.0) then
c         ii=ii+1
c         krelm(karis,ii)=kg(icopo,jcopo)
c        end if
c        if(kg(jcopo,icopo).ne.0) then
c         ii=ii+1
c         krelm(karis,ii)=kg(jcopo,icopo)
c        end if
c       end if
c      end do
c     end do
C
      karis=0                          ! new form
      do icoel=1,ncoel
       karis=karis+1
       krnoe(karis,1)=lnodc(1,icoel)
       krnoe(karis,2)=lnodc(2,icoel)
       krelm(karis,1)=icoel
       krelm(karis,2)=0
       karis=karis+1
       krnoe(karis,1)=lnodc(2,icoel)
       krnoe(karis,2)=lnodc(3,icoel)
       krelm(karis,1)=icoel
       krelm(karis,2)=0
       karis=karis+1
       krnoe(karis,1)=lnodc(3,icoel)
       krnoe(karis,2)=lnodc(1,icoel)
       krelm(karis,1)=icoel
       krelm(karis,2)=0
      end do
c
      naris=karis
      do iaris=karis-1,1,-1
       do jaris=iaris+1,naris
        if((krnoe(jaris,1).eq.krnoe(iaris,2)).and.
     .     (krnoe(jaris,2).eq.krnoe(iaris,1))) then
         krelm(iaris,2)=krelm(jaris,1)
         do laris=jaris+1,naris
          krnoe(laris-1,1)=krnoe(laris,1)
          krnoe(laris-1,2)=krnoe(laris,2)
          krelm(laris-1,1)=krelm(laris,1)
          krelm(laris-1,2)=krelm(laris,2)
         end do
         naris=naris-1
        end if
       end do
      end do
      karis=naris
C
C**** REDEFINES KRNOE ACCORDING TO CONNECTIVITY OF ELEMENT KRELM(*,1)
C     (THIS IS NECESSARY FOR THE DIAGONAL INVERSION OPERATION)
C
      do iaris=1,karis
       if((krnoe(iaris,1).eq.lnodc(3,krelm(iaris,1)).and.
     .     krnoe(iaris,2).eq.lnodc(2,krelm(iaris,1))).or.
     .    (krnoe(iaris,1).eq.lnodc(2,krelm(iaris,1)).and.
     .     krnoe(iaris,2).eq.lnodc(1,krelm(iaris,1))).or.
     .    (krnoe(iaris,1).eq.lnodc(1,krelm(iaris,1)).and.
     .     krnoe(iaris,2).eq.lnodc(3,krelm(iaris,1)))) then
        ii=krnoe(iaris,1)
        krnoe(iaris,1)=krnoe(iaris,2)
        krnoe(iaris,2)=ii
       end if
      end do
C
      if(karis.gt.maris) then
c      CALL RUNENDF('ERROR: NUMBER OF EDGES GT MAXIMUM')
       write(83,*) 'karis gt maris - karis=',karis           ! mlit-mesh
       stop
      end if
c
c**** Computes:
c
c     kelar: edge numbers for each contour element (respecting the
c            original connectivity)
c
      do icoel=1,ncoel
       do iaris=1,karis
        if(((krnoe(iaris,1).eq.lnodc(1,icoel).and.
     .       krnoe(iaris,2).eq.lnodc(2,icoel))).or.
     .     ((krnoe(iaris,1).eq.lnodc(2,icoel).and.
     .       krnoe(iaris,2).eq.lnodc(1,icoel)))) kelar(icoel,1)=iaris
        if(((krnoe(iaris,1).eq.lnodc(2,icoel).and.
     .       krnoe(iaris,2).eq.lnodc(3,icoel))).or.
     .     ((krnoe(iaris,1).eq.lnodc(3,icoel).and.
     .       krnoe(iaris,2).eq.lnodc(2,icoel)))) kelar(icoel,2)=iaris
        if(((krnoe(iaris,1).eq.lnodc(3,icoel).and.
     .       krnoe(iaris,2).eq.lnodc(1,icoel))).or.
     .     ((krnoe(iaris,1).eq.lnodc(1,icoel).and.
     .       krnoe(iaris,2).eq.lnodc(3,icoel)))) kelar(icoel,3)=iaris
       end do
      end do
C
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE DIAGON(LNODN,NNUBE,NCPNU,LELNU,KRELM,KRNOE,KELAR,
     .                  XNN,COSAM,
     .                  MAXCE,MAXCP,MELCP,MARIS,
     .                  NDIMEF,NNCEL,KARIS,LURESF)
C***********************************************************************
C
C**** THIS ROUTINE IMPROVES CONTOUR ELEMENT QUALITY BY INVERTING ITS
C     DIAGONAL ACCORDING TO THE FOLLOWING CRITERIA:
C     - ANGLE
C     - QUALITY
C     - CONVEXITY
C     - AREA RATIO
C
C***********************************************************************
      implicit double precision (a-h,o-z)
C
      dimension lnodn(nncel+1,maxce)
      dimension nnube(maxcp,melcp*nncel),                  ! work arrays
     .          ncpnu(maxcp,5),
     .          lelnu(maxcp,melcp),
     .          krelm(maris,2), krnoe(maris,2),
     .          kelar(maxce,3)
      dimension xnn(maxcp,2*ndimef)
      DIMENSION edge1(3), edge2(3), edge3(3), areai(3), areaj(3)
C
      icrit=2            ! two quality criteria (better as input)
      do iaris=1,karis
       ie1=-1                          ! to deal with non-conforming el
       ie2=-1
       if(krelm(iaris,2).ne.0) then
        ie1=lnodn(nncel+1,krelm(iaris,1))
        ie2=lnodn(nncel+1,krelm(iaris,2))
       end if
       if((ie1.eq.0.and.ie2.eq.0).or.  ! only (i & b) internal edges
     .    (ie1.eq.1.and.ie2.eq.1)) then
        i1=lnodn(1,krelm(iaris,1))     ! first element sharing the edge
        i2=lnodn(2,krelm(iaris,1))
        i3=lnodn(3,krelm(iaris,1))
        do j=1,3                       ! 3D
         edge1(j)=xnn(i2,j)-xnn(i1,j)
         edge2(j)=xnn(i3,j)-xnn(i1,j)
         edge3(j)=xnn(i3,j)-xnn(i2,j)
        end do
        call vecpro(3,edge1,edge2,areai)
        areaim=areai(1)*areai(1)+areai(2)*areai(2)+areai(3)*areai(3)
        do i=1,3                       ! normalization
         areai(i)=areai(i)/dsqrt(areaim)
        end do
        areaim=0.5d0*dsqrt(areaim)
        suarea=areaim
        edgem1=edge1(1)*edge1(1)+edge1(2)*edge1(2)+edge1(3)*edge1(3)
        edgem1=dsqrt(edgem1)
        edgem2=edge2(1)*edge2(1)+edge2(2)*edge2(2)+edge2(3)*edge2(3)
        edgem2=dsqrt(edgem2)
        edgem3=edge3(1)*edge3(1)+edge3(2)*edge3(2)+edge3(3)*edge3(3)
        edgem3=dsqrt(edgem3)
        if(icrit.eq.2) then
         edmax=edgem1
         if(edgem2.gt.edmax) edmax=edgem2
         if(edgem3.gt.edmax) edmax=edgem3
         rinsc=2.0d0*areaim/(edgem1+edgem2+edgem3)
         qual1=rinsc/edmax
         quali=qual1
        end if
c
        j1=lnodn(1,krelm(iaris,2))     ! second element sharing the edge
        j2=lnodn(2,krelm(iaris,2))
        j3=lnodn(3,krelm(iaris,2))
        do j=1,3                       ! 3D
         edge1(j)=xnn(j2,j)-xnn(j1,j)
         edge2(j)=xnn(j3,j)-xnn(j1,j)
         edge3(j)=xnn(j3,j)-xnn(j2,j)
        end do
        call vecpro(3,edge1,edge2,areaj)
        areajm=areaj(1)*areaj(1)+areaj(2)*areaj(2)+areaj(3)*areaj(3)
        do i=1,3                       ! normalization
         areaj(i)=areaj(i)/dsqrt(areajm)
        end do
        areajm=0.5d0*dsqrt(areajm)
        suarea=suarea+areajm
        edgem1=edge1(1)*edge1(1)+edge1(2)*edge1(2)+edge1(3)*edge1(3)
        edgem1=dsqrt(edgem1)
        edgem2=edge2(1)*edge2(1)+edge2(2)*edge2(2)+edge2(3)*edge2(3)
        edgem2=dsqrt(edgem2)
        edgem3=edge3(1)*edge3(1)+edge3(2)*edge3(2)+edge3(3)*edge3(3)
        edgem3=dsqrt(edgem3)
        if(icrit.eq.2) then
         edmax=edgem1
         if(edgem2.gt.edmax) edmax=edgem2
         if(edgem3.gt.edmax) edmax=edgem3
         rinsc=2.0d0*areajm/(edgem1+edgem2+edgem3)
         qual2=rinsc/edmax
         if(qual2.lt.quali) quali=qual2
        end if
        aratio=areaim/areajm                  ! checks area ratio
        if(aratio.gt.1.0d0) aratio=1.0d0/aratio
c
        cosarea=areai(1)*areaj(1)+areai(2)*areaj(2)+areai(3)*areaj(3)
        if(cosarea.gt.cosam) then             ! nearly coplanar elements
c
         k1=krelm(iaris,1)
         k2=krelm(iaris,2)
c
         kk=0
         if((lnodn(2,k1).eq.lnodn(1,k2)).and.
     .      (lnodn(3,k1).eq.lnodn(3,k2))) then  ! abc - bdc -> adc - abd
          ia1=lnodn(1,k1)                       ! a
          ia2=lnodn(2,k1)                       ! b
          ia3=lnodn(3,k1)                       ! c
          ja3=lnodn(2,k2)                       ! d
          kk=kk+1
         end if
         if((lnodn(2,k1).eq.lnodn(3,k2)).and.
     .      (lnodn(3,k1).eq.lnodn(2,k2))) then  ! abc - dcb -> adc - abd
          ia1=lnodn(1,k1)                       ! a
          ia2=lnodn(2,k1)                       ! b
          ia3=lnodn(3,k1)                       ! c
          ja3=lnodn(1,k2)                       ! d
          kk=kk+1
         end if
         if((lnodn(2,k1).eq.lnodn(2,k2)).and.
     .      (lnodn(3,k1).eq.lnodn(1,k2))) then  ! abc - cbd -> adc - abd
          ia1=lnodn(1,k1)                       ! a
          ia2=lnodn(2,k1)                       ! b
          ia3=lnodn(3,k1)                       ! c
          ja3=lnodn(3,k2)                       ! d
          kk=kk+1
         end if
         if((lnodn(1,k1).eq.lnodn(1,k2)).and.
     .      (lnodn(2,k1).eq.lnodn(3,k2))) then  ! bca - bdc -> adc - abd
          ia1=lnodn(3,k1)                       ! a
          ia2=lnodn(1,k1)                       ! b
          ia3=lnodn(2,k1)                       ! c
          ja3=lnodn(2,k2)                       ! d
          kk=kk+1
         end if
         if((lnodn(1,k1).eq.lnodn(3,k2)).and.
     .      (lnodn(2,k1).eq.lnodn(2,k2))) then  ! bca - dcb -> adc - abd
          ia1=lnodn(3,k1)                       ! a
          ia2=lnodn(1,k1)                       ! b
          ia3=lnodn(2,k1)                       ! c
          ja3=lnodn(1,k2)                       ! d
          kk=kk+1
         end if
         if((lnodn(1,k1).eq.lnodn(2,k2)).and.
     .      (lnodn(2,k1).eq.lnodn(1,k2))) then  ! bca - cbd -> adc - abd
          ia1=lnodn(3,k1)                       ! a
          ia2=lnodn(1,k1)                       ! b
          ia3=lnodn(2,k1)                       ! c
          ja3=lnodn(3,k2)                       ! d
          kk=kk+1
         end if
         if((lnodn(3,k1).eq.lnodn(1,k2)).and.
     .      (lnodn(1,k1).eq.lnodn(3,k2))) then  ! cab - bdc -> adc - abd
          ia1=lnodn(2,k1)                       ! a
          ia2=lnodn(3,k1)                       ! b
          ia3=lnodn(1,k1)                       ! c
          ja3=lnodn(2,k2)                       ! d
          kk=kk+1
         end if
         if((lnodn(3,k1).eq.lnodn(3,k2)).and.
     .      (lnodn(1,k1).eq.lnodn(2,k2))) then  ! cab - dcb -> adc - abd
          ia1=lnodn(2,k1)                       ! a
          ia2=lnodn(3,k1)                       ! b
          ia3=lnodn(1,k1)                       ! c
          ja3=lnodn(1,k2)                       ! d
          kk=kk+1
         end if
         if((lnodn(3,k1).eq.lnodn(2,k2)).and.
     .      (lnodn(1,k1).eq.lnodn(1,k2))) then  ! cab - cbd -> adc - abd
          ia1=lnodn(2,k1)                       ! a
          ia2=lnodn(3,k1)                       ! b
          ia3=lnodn(1,k1)                       ! c
          ja3=lnodn(3,k2)                       ! d
          kk=kk+1
         end if
         if(kk.ne.1) then
          write(luresf,*) 'k1,k2,kk=',k1,k2,kk
          write(luresf,*) 'i1,i2,i3=',i1,i2,i3
          write(luresf,*) 'j1,j2,j3=',j1,j2,j3
c         call runendf('error: kk ne 1 in mlirt3d.f')
          write(83,*) 'error: kk ne 1 in mlirt3d.f'          ! mlit-mesh
          stop
         end if
c
         do j=1,3                      ! 3D
          edge1(j)=xnn(ja3,j)-xnn(ia1,j)
          edge2(j)=xnn(ia3,j)-xnn(ia1,j)
          edge3(j)=xnn(ia3,j)-xnn(ja3,j)
         end do
         call vecpro(3,edge1,edge2,areai)
         areaim=areai(1)*areai(1)+areai(2)*areai(2)+areai(3)*areai(3)
         areaim=0.5d0*dsqrt(areaim)
         suareaa=areaim
         edgem1=edge1(1)*edge1(1)+edge1(2)*edge1(2)+edge1(3)*edge1(3)
         edgem1=dsqrt(edgem1)
         edgem2=edge2(1)*edge2(1)+edge2(2)*edge2(2)+edge2(3)*edge2(3)
         edgem2=dsqrt(edgem2)
         edgem3=edge3(1)*edge3(1)+edge3(2)*edge3(2)+edge3(3)*edge3(3)
         edgem3=dsqrt(edgem3)
         if(icrit.eq.2) then
          edmax=edgem1
          if(edgem2.gt.edmax) edmax=edgem2
          if(edgem3.gt.edmax) edmax=edgem3
          rinsc=2.0d0*areaim/(edgem1+edgem2+edgem3)
          qual1=rinsc/edmax
          qualia=qual1
         end if
c
         do j=1,3                      ! 3D
          edge1(j)=xnn(ia2,j)-xnn(ia1,j)
          edge2(j)=xnn(ja3,j)-xnn(ia1,j)
          edge3(j)=xnn(ja3,j)-xnn(ia2,j)
         end do
         call vecpro(3,edge1,edge2,areaj)
         areajm=areaj(1)*areaj(1)+areaj(2)*areaj(2)+areaj(3)*areaj(3)
         areajm=0.5d0*dsqrt(areajm)
         suareaa=suareaa+areajm
         edgem1=edge1(1)*edge1(1)+edge1(2)*edge1(2)+edge1(3)*edge1(3)
         edgem1=dsqrt(edgem1)
         edgem2=edge2(1)*edge2(1)+edge2(2)*edge2(2)+edge2(3)*edge2(3)
         edgem2=dsqrt(edgem2)
         edgem3=edge3(1)*edge3(1)+edge3(2)*edge3(2)+edge3(3)*edge3(3)
         edgem3=dsqrt(edgem3)
         if(icrit.eq.2) then
          edmax=edgem1
          if(edgem2.gt.edmax) edmax=edgem2
          if(edgem3.gt.edmax) edmax=edgem3
          rinsc=2.0d0*areajm/(edgem1+edgem2+edgem3)
          qual2=rinsc/edmax
          if(qual2.lt.qualia) qualia=qual2
         end if
         aratioa=areaim/areajm                   ! checks area ratio
         if(aratioa.gt.1.0d0) aratioa=1.0d0/aratioa
C
         psuare=dabs((suarea-suareaa)/suarea)    ! checks convexity
C
         if(qualia.gt.quali*1.001d0.and.         ! do diagonal inversion
     .      psuare.lt.1.0d-03.and.aratioa.gt.0.5d0*aratio) then
c
          lnodn(1,krelm(iaris,1))=ia1  ! a
          lnodn(2,krelm(iaris,1))=ja3  ! d
          lnodn(3,krelm(iaris,1))=ia3  ! c
          lnodn(1,krelm(iaris,2))=ia1  ! a
          lnodn(2,krelm(iaris,2))=ia2  ! b
          lnodn(3,krelm(iaris,2))=ja3  ! d
C
          krnoe(iaris,1)=ia1           ! a
          krnoe(iaris,2)=ja3           ! d
c
c**** Redefines kelar
c
          do laris=1,3
           jaris=kelar(krelm(iaris,2),laris) ! old kelar
           if((krnoe(jaris,1).eq.ja3.and.    ! d & c
     .         krnoe(jaris,2).eq.ia3)) krelm(jaris,1)=krelm(iaris,1)
           if((krnoe(jaris,1).eq.ia3.and.    ! c & d
     .         krnoe(jaris,2).eq.ja3)) krelm(jaris,2)=krelm(iaris,1)
c
           if(((krnoe(jaris,1).eq.ja3.and.   ! d & c
     .          krnoe(jaris,2).eq.ia3)).or.
     .        ((krnoe(jaris,1).eq.ia3.and.   ! c & d
     .          krnoe(jaris,2).eq.ja3))) i1ar2=jaris
           if(((krnoe(jaris,1).eq.ia2.and.   ! b & d
     .          krnoe(jaris,2).eq.ja3)).or.
     .        ((krnoe(jaris,1).eq.ja3.and.   ! d & b
     .          krnoe(jaris,2).eq.ia2))) i2ar3=jaris
c
           jaris=kelar(krelm(iaris,1),laris) ! old kelar
           if((krnoe(jaris,1).eq.ia1.and.    ! a & b
     .         krnoe(jaris,2).eq.ia2)) krelm(jaris,1)=krelm(iaris,2)
           if((krnoe(jaris,1).eq.ia2.and.    ! b & a
     .         krnoe(jaris,2).eq.ia1)) krelm(jaris,2)=krelm(iaris,2)
c
           if(((krnoe(jaris,1).eq.ia1.and.   ! a & b
     .          krnoe(jaris,2).eq.ia2)).or.
     .        ((krnoe(jaris,1).eq.ia2.and.   ! b & a
     .          krnoe(jaris,2).eq.ia1))) i2ar2=jaris
           if(((krnoe(jaris,1).eq.ia3.and.   ! c & a
     .          krnoe(jaris,2).eq.ia1)).or.
     .        ((krnoe(jaris,1).eq.ia1.and.   ! a & c
     .          krnoe(jaris,2).eq.ia3))) i1ar3=jaris
          end do
c
          kelar(krelm(iaris,1),1)=iaris      ! redefines kelar
          kelar(krelm(iaris,1),2)=i1ar2
          kelar(krelm(iaris,1),3)=i1ar3
          kelar(krelm(iaris,2),1)=iaris
          kelar(krelm(iaris,2),2)=i2ar2
          kelar(krelm(iaris,2),3)=i2ar3
c
c**** Redefines ncpnu, nnube & lelnu
c
          ncpnu(ia1,1)=ncpnu(ia1,1)+1        ! a
          ncpnu(ia1,2)=ncpnu(ia1,2)+1
          nnube(ia1,ncpnu(ia1,1))=ja3 
          lelnu(ia1,ncpnu(ia1,2))=krelm(iaris,2)
c
          do kn=1,ncpnu(ia2,1)               ! b
           if(nnube(ia2,kn).eq.ia3) nnube(ia2,kn)=0
          end do
          jn=0
          do kn=1,ncpnu(ia2,1)
           if(nnube(ia2,kn).ne.0) then
            jn=jn+1
            nnube(ia2,jn)=nnube(ia2,kn)
           end if
          end do
          do kn=1,ncpnu(ia2,2)
           if(lelnu(ia2,kn).eq.krelm(iaris,1)) lelnu(ia2,kn)=0
          end do
          jn=0
          do kn=1,ncpnu(ia2,2)
           if(lelnu(ia2,kn).ne.0) then
            jn=jn+1
            lelnu(ia2,jn)=lelnu(ia2,kn)
           end if
          end do
          ncpnu(ia2,1)=ncpnu(ia2,1)-1
          ncpnu(ia2,2)=ncpnu(ia2,2)-1
c
          do kn=1,ncpnu(ia3,1)               ! c
           if(nnube(ia3,kn).eq.ia2) nnube(ia3,kn)=0
          end do
          jn=0
          do kn=1,ncpnu(ia3,1)
           if(nnube(ia3,kn).ne.0) then
            jn=jn+1
            nnube(ia3,jn)=nnube(ia3,kn)
           end if
          end do
          do kn=1,ncpnu(ia3,2)
           if(lelnu(ia3,kn).eq.krelm(iaris,2)) lelnu(ia3,kn)=0
          end do
          jn=0
          do kn=1,ncpnu(ia3,2)
           if(lelnu(ia3,kn).ne.0) then
            jn=jn+1
            lelnu(ia3,jn)=lelnu(ia3,kn)
           end if
          end do
          ncpnu(ia3,1)=ncpnu(ia3,1)-1
          ncpnu(ia3,2)=ncpnu(ia3,2)-1
c
          ncpnu(ja3,1)=ncpnu(ja3,1)+1        ! d
          ncpnu(ja3,2)=ncpnu(ja3,2)+1
          nnube(ja3,ncpnu(ja3,1))=ia1 
          lelnu(ja3,ncpnu(ja3,2))=krelm(iaris,1)
c
         end if          ! qualia gt quali & psuare gt 1.0d-03
        end if           ! cosarea gt cosam
       end if            ! ie1=0 & ie2=0 or ...
      end do             ! iaris=1,karis
c
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE OUTGID(LNODN,XNN,CCONT,
     .                  MAXCE,NNCEL,NDIMEF,KMOVI,KPREC,MAXCP,
     .                  NCOEN,NCOPN,ISTEPF,INDEX)
C***********************************************************************
C
C**** THIS ROUTINE PRINTS THE CONTOUR MESH FOR GiD & RESTART
C
C***********************************************************************
      implicit double precision (a-h,o-z)
C
      dimension lnodn(nncel+1,maxce)
      dimension xnn(maxcp,*)              ! to allow the use of xn & xnn
      DIMENSION CCONT(NDIMEF+2+KMOVI*(NDIMEF+KPREC),MAXCP)
c
c***********************************************************************
c     GiD for output mesh
c
      write(83,*) 'titulo: problema - paso=',ISTEPF
      write(83,*) 'subtitulo: problema'
      write(83,*) 'a'
      write(83,*) 'a'
      write(83,*) 'a'
c
      write(83,*) 'c'
      nnodx=3
      write(83,*) ncoen*3, ncopn, nnodx
c
      write(83,*) 'c'
      do ipoin=1,ncopn
       if(index.eq.0) then            ! call from conmov.f
        x=ccont(1,ipoin)
        y=ccont(2,ipoin)
        z=ccont(3,ipoin)
       else                           ! call from mlirt3d.f
        x=xnn(ipoin,1)
        y=xnn(ipoin,2)
        z=xnn(ipoin,3)
       end if
       write(83,901) ipoin,x,y,z      ! formatted
      end do
      write(83,*) 'c'
c
      do ielem=1,ncoen
       write(83,902) ielem,(lnodn(i,ielem),i=1,3),lnodn(1,ielem),
     .               lnodn(4,ielem)+1
       write(83,902) ielem+ncoen,(lnodn(i,ielem),i=1,3),lnodn(2,ielem),
     .               lnodn(4,ielem)+1
       write(83,902) ielem+2*ncoen,(lnodn(i,ielem),i=1,3),
     .                                                  lnodn(3,ielem),
     .               lnodn(4,ielem)+1
      end do
C
c***********************************************************************
c     Restart for output mesh
c
      write(83,*) ' '
      write(83,*) 'CONTOUR'
      write(83,*) 'COORDINATES'
      do ipoin=1,ncopn
       if(index.eq.0) then            ! call from conmov.f
        x=ccont(1,ipoin)
        y=ccont(2,ipoin)
        z=ccont(3,ipoin)
       else                           ! call from mlirt3d.f
        x=xnn(ipoin,1)
        y=xnn(ipoin,2)
        z=xnn(ipoin,3)
       end if
       ie=CCONT(NDIMEF+2,IPOIN)       ! reference element to search
       ix=CCONT(NDIMEF+3,IPOIN)       ! prescribed x velocity
       iy=CCONT(NDIMEF+4,IPOIN)       ! prescribed y velocity
       iz=CCONT(NDIMEF+5,IPOIN)       ! prescribed z velocity
       ip=CCONT(NDIMEF+6,IPOIN)       ! prescribed pressure
       write(83,903) ipoin,x,y,z,ie,ix,iy,iz,ip     ! formatted
      end do
      write(83,*) 'END_COORDINATES'

      write(83,*) 'ELEMENTS'
      do ielem=1,ncoen
       write(83,902) ielem,(lnodn(i,ielem),i=1,3),lnodn(4,ielem)
      end do
      write(83,*) 'END_ELEMENTS'
      write(83,*) 'END_CONTOUR'
      write(83,*) ' '
C
  901 format(i8,3e15.6)
  902 format(i8,4i8,i5)
  903 format(i8,3e15.6,i8,4i3)
C
c***********************************************************************
c
c     write(luresf,*) 'COORDINATES'
c     do ipoin=1,ncopn
c      x=xnn(ipoin,1)
c      y=xnn(ipoin,2)
c      z=xnn(ipoin,3)
c      write(luresf,*) ipoin,x,y,z,'  1.0   1.0   1.0   1.0   0.0'
c     end do
c     write(luresf,*) 'END_COORDINATES'
c
c     write(luresf,*)'ELEMENTS'
c     do ielem=1,ncoen
c      write(luresf,902) ielem,(lnodn(i,ielem),i=1,3)
c     end do
c     write(luresf,*)'END_ELEMENTS'
c
c***********************************************************************
C
      RETURN
      END
C
c***********************************************************************
      SUBROUTINE REPOSB(LNODN,NNUBE,NCPNU,LELNU,KRELM,KRNOE,
     .                  XNN,COSAM,WEIML,
     .                  NCOPN,
     .                  MAXCE,MAXCP,MELCP,MARIS,
     .                  NDIMEF,NNCEL,KARIS,LURESF)
C***********************************************************************C
C**** THIS ROUTINE PERFORMS PATCH-BASED & EDGE-BASED BARICENTRIC
C     REPOSITIONING
C
C***********************************************************************
      implicit double precision (a-h,o-z)
C
      dimension lnodn(nncel+1,maxce)
      dimension nnube(maxcp,melcp*nncel),                  ! work arrays
     .          ncpnu(maxcp,5),
     .          lelnu(maxcp,melcp),
     .          krelm(maris,2), krnoe(maris,2)
      dimension xnn(maxcp,2*ndimef)
      DIMENSION edge1(3), edge2(3), edge3(3), areai(3), areaj(3),
     .          xpnew(3)
C
      do icopo=1,ncopn
       if(ncpnu(icopo,5).eq.0) then     ! only interface internal cp
        nelcp=ncpnu(icopo,2)            ! computes normals
        ku=lelnu(icopo,1)               ! first element of the cloud
        i1=lnodn(1,ku)                  ! normal to this element
        i2=lnodn(2,ku)
        i3=lnodn(3,ku)
        do j=1,3                        ! 3D
         edge1(j)=xnn(i2,j)-xnn(i1,j)
         edge2(j)=xnn(i3,j)-xnn(i1,j)
        end do
        call vecpro(3,edge1,edge2,areai)
        areaim=areai(1)*areai(1)+areai(2)*areai(2)+areai(3)*areai(3)
        do i=1,3                        ! normalization
         areai(i)=areai(i)/dsqrt(areaim)
        end do
        kpr=1
        if(nelcp.gt.1) then
         do ielcp=2,nelcp               ! normal to the other elements
          ku=lelnu(icopo,ielcp)
          i1=lnodn(1,ku)
          i2=lnodn(2,ku)
          i3=lnodn(3,ku)
          do j=1,3                      ! 3D
           edge1(j)=xnn(i2,j)-xnn(i1,j)
           edge2(j)=xnn(i3,j)-xnn(i1,j)
          end do
          call vecpro(3,edge1,edge2,areaj)
          areajm=areaj(1)*areaj(1)+areaj(2)*areaj(2)+areaj(3)*areaj(3)
          do i=1,3                      ! normalization
           areaj(i)=areaj(i)/dsqrt(areajm)
          end do
          cosarea=areai(1)*areaj(1)+areai(2)*areaj(2)+areai(3)*areaj(3)
          if(cosarea.gt.cosam) kpr=kpr+1      ! nearly coplanar elements
         end do   ! ielcp
        end if
        weix=weiml                      ! to avoid excessive smoothing
        if(kpr.eq.nelcp) weix=0.5d0     ! default value
c
        do idime=1,ndimef
         xpnew(idime)=0.0d0
        end do
        do ix=1,ncpnu(icopo,1)
         if(icopo.ne.nnube(icopo,ix)) then
          do idime=1,ndimef
           xpnew(idime)=xpnew(idime)+xnn(nnube(icopo,ix),idime)
          end do
         end if
        end do
        do idime=1,ndimef               ! apply baricentric reposit.
         xnn(icopo,idime)=(1.0d0-weix)*xnn(icopo,idime)+
     .                              weix*xpnew(idime)/(ncpnu(icopo,1)-1)
        end do
       end if    ! ncpnu(icopo,5)=0
C
       if(ncpnu(icopo,5).eq.2) then     ! only boundary internal cp
        nelcp=ncpnu(icopo,2)            ! computes normals
        ku=lelnu(icopo,1)               ! first element of the cloud
        i1=lnodn(1,ku)                  ! normal to this element
        i2=lnodn(2,ku)
        i3=lnodn(3,ku)
        do j=1,3                        ! 3D
         edge1(j)=xnn(i2,j)-xnn(i1,j)
         edge2(j)=xnn(i3,j)-xnn(i1,j)
        end do
        call vecpro(3,edge1,edge2,areai)
        areaim=areai(1)*areai(1)+areai(2)*areai(2)+areai(3)*areai(3)
        do i=1,3                        ! normalization
         areai(i)=areai(i)/dsqrt(areaim)
        end do
        kpr=1
        if(nelcp.gt.1) then
         do ielcp=2,nelcp               ! normal to the other elements
          ku=lelnu(icopo,ielcp)
          i1=lnodn(1,ku)
          i2=lnodn(2,ku)
          i3=lnodn(3,ku)
          do j=1,3                      ! 3D
           edge1(j)=xnn(i2,j)-xnn(i1,j)
           edge2(j)=xnn(i3,j)-xnn(i1,j)
          end do
          call vecpro(3,edge1,edge2,areaj)
          areajm=areaj(1)*areaj(1)+areaj(2)*areaj(2)+areaj(3)*areaj(3)
          do i=1,3                      ! normalization
           areaj(i)=areaj(i)/dsqrt(areajm)
          end do
          cosarea=areai(1)*areaj(1)+areai(2)*areaj(2)+areai(3)*areaj(3)
          if(cosarea.gt.cosam) kpr=kpr+1      ! nearly coplanar elements
         end do   ! ielcp
        end if
        weix=weiml                      ! to avoid excessive smoothing
        if(kpr.eq.nelcp) weix=0.5d0     ! default value
c
        do idime=1,ndimef
         xpnew(idime)=0.0d0
        end do
        do ix=1,ncpnu(icopo,1)
         if(icopo.ne.nnube(icopo,ix)) then
          do idime=1,ndimef
           xpnew(idime)=xpnew(idime)+xnn(nnube(icopo,ix),idime)
          end do
         end if
        end do
        do idime=1,ndimef               ! apply baricentric reposit.
         if(int(xnn(icopo,idime+ndimef)).eq.0)
     .    xnn(icopo,idime)=(1.0d0-weix)*xnn(icopo,idime)+
     .                              weix*xpnew(idime)/(ncpnu(icopo,1)-1)
        end do
       end if    ! ncpnu(icopo,5)=2
      end do     ! icopo
C
C**** STEP 8: EDGE-BASED BARICENTRIC REPOSITIONING
C
C     Note: for boundary contour points, krelm(iaris,2) ne 0
C
      do icopo=1,ncopn
       if(ncpnu(icopo,5).eq.1) then     ! only boundary contour points
        do icop1=1,ncpnu(icopo,1)-1
         jcopo=nnube(icopo,icop1)
         if(icopo.ne.jcopo.and.ncpnu(jcopo,5).eq.1) then
          do icop2=icop1+1,ncpnu(icopo,1)
           kcopo=nnube(icopo,icop2)
           if(icopo.ne.kcopo.and.ncpnu(kcopo,5).eq.1) then
            i1=icopo
            i2=jcopo
            i3=kcopo
            ix1=0                       ! check bcp belong to b edges
            ix2=0
            do iaris=1,karis
             if((i1.eq.krnoe(iaris,1).and.i2.eq.krnoe(iaris,2)).or.
     .          (i1.eq.krnoe(iaris,2).and.i2.eq.krnoe(iaris,1))) then
              ie1=lnodn(nncel+1,krelm(iaris,1))
              ie2=lnodn(nncel+1,krelm(iaris,2))
              if((ie1.eq.1.and.ie2.eq.0).or.
     .           (ie1.eq.0.and.ie2.eq.1)) ix1=1
             end if
             if((i1.eq.krnoe(iaris,1).and.i3.eq.krnoe(iaris,2)).or.
     .          (i1.eq.krnoe(iaris,2).and.i3.eq.krnoe(iaris,1))) then
              ie1=lnodn(nncel+1,krelm(iaris,1))
              ie2=lnodn(nncel+1,krelm(iaris,2))
              if((ie1.eq.1.and.ie2.eq.0).or.
     .           (ie1.eq.0.and.ie2.eq.1)) ix2=1
             end if
            end do
            if(ix1.eq.1.and.
     .         ix2.eq.1) then           ! only bcp belonging to b edges
             do j=1,3                   ! computes normals
              edge1(j)=xnn(i2,j)-xnn(i1,j)
              edge2(j)=xnn(i3,j)-xnn(i1,j)
             end do
             edge1l=edge1(1)*edge1(1)+edge1(2)*edge1(2)+
     .              edge1(3)*edge1(3)
             edge2l=edge2(1)*edge2(1)+edge2(2)*edge2(2)+
     .              edge2(3)*edge2(3)
             do i=1,3                   ! normalization
              edge1(i)=edge1(i)/dsqrt(edge1l)
              edge2(i)=edge2(i)/dsqrt(edge2l)
             end do
             cosedg=edge1(1)*edge2(1)+edge1(2)*edge2(2)+
     .              edge1(3)*edge2(3)
             cosedg=-cosedg             ! angles > 90 (close to 180)
             weix=weiml                 ! to avoid excessive smoothing
             if(cosedg.gt.cosam) weix=0.5d0            ! default value
c
             do idime=1,ndimef          ! apply baricentric reposit.
              if(int(xnn(icopo,idime+ndimef)).eq.0)
     .         xnn(icopo,idime)=(1.0d0-weix)*xnn(icopo,idime)+
     .                    weix*(xnn(jcopo,idime)+xnn(kcopo,idime))/2.0d0
             end do
            end if
           end if
          end do
         end if
        end do
       end if    ! ncpnu(icopo,5)=1
      end do     ! icopo
C
      RETURN
      END
C
C***********************************************************************
      SUBROUTINE VECPRO(N,V1,V2,V3)
C***********************************************************************
C
C**** TRIDIMENSIONAL VECTORIAL PRODUCT OF TWO VECTORS  V3 = V1 X V2
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION V1(N),V2(N),V3(N)
C
      V3(1) = V1(2)*V2(3) - V1(3)*V2(2)
      V3(2) = V1(3)*V2(1) - V1(1)*V2(3)
      V3(3) = V1(1)*V2(2) - V1(2)*V2(1)
      RETURN
      END
C***********************************************************************
      SUBROUTINE INVERT(A,NMAX,NDM)
C***********************************************************************
C
C**** THIS ROUTINE PERFORMS THE INVERSION OF A NDM*NDM SQUARE MATRIX 
C     OR JUST PART OF IT (NMAX*NMAX)
C
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NDM,NDM)
      DO 200 N = 1,NMAX
      D = A(N,N)
      DO 100 J = 1,NMAX
100   A(N,J) = -A(N,J)/D
      DO 150 I = 1,NMAX
      IF(N.EQ.I) GO TO 150
      DO 140 J = 1,NMAX
      IF(N.NE.J) A(I,J) = A(I,J) +A(I,N)*A(N,J)
140   CONTINUE
150   A(I,N) = A(I,N)/D
      A(N,N) = 1.0/D
200   CONTINUE
      RETURN
      END
C***********************************************************************
