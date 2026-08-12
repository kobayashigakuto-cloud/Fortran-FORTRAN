C-------declaration of variables--------------
C
        implicit none
        INTEGER TNST,NST,stepN,timeinit,DST
        INTEGER O,H
        INTEGER I,J,K,Q
        REAL(8) dx,dy,dz,x,y,z,az,adz
        REAL(8) dr,dist,dlimit,datoms
        INTEGER Nstart,dstep
        INTEGER NumC,num
        INTEGER divx,divy,divz
        REAL(8) space     
        CHARACTER(4) ATOM(10)
        INTEGER NTION,NCOMPO,STP,NION(10),IONS(2,10)
        INTEGER IPV(3,5000),preIPV(3,5000)
        REAL(8) CELL(3),P(3,2160,50000),dcell(3,50000)
C-------open files
C
        OPEN(11,FILE='file07.dat',STATUS='old')
        OPEN(12,file='file09p.dat',status='old')
        OPEN(22,file='spaceH2O.dat',status='replace')
        OPEN(31,FILE='parardf.dat',STATUS='old')
C
C-------reading of 'file07'
C
        READ (11,101) NTION,NCOMPO
        READ (11,102) (ATOM(I),I=1,NCOMPO)
        READ (11,103) (NION(I),I=1,NCOMPO)
        READ (11,103) (IONS(1,I),I=1,NCOMPO)
        READ (11,103) (IONS(2,I),I=1,NCOMPO)
  101   FORMAT (/I7,I3)
  102   FORMAT (10(2X,A4))
  103   FORMAT (10I6)
C
c-------reading of 'parardf.dat'
C
        READ (31,104) timeinit,dr,stepN,STP
  104   FORMAT(I8/F8.2/I8/I8)
C          
C-------reading of 'file09p'
C
      DO 210 I=1,STP
        READ (12,105) (CELL(K),K=1,3)
        READ (12,106) ((IPV(K,J),K=1,3),J=1,NTION)
  105   FORMAT (10X,F7.3,21X,F7.3,21X,F7.3)
  106   FORMAT (18I4)
            dcell(1,I)=CELL(1)
            dcell(2,I)=CELL(2)
            dcell(3,I)=CELL(3)
           DO 650 J=1,NTION
           DO 660 K=1,3
               P(K,J,I)=CELL(K)*dble(IPV(K,J))/9d3
  660      END DO
  650     END DO
  210     END DO
C
C-------rdf(O atom around O atom)
C
      az=0
      DO I=timeinit+9,STP,10
       az=az+dcell(3,I)-dcell(3,STP)
      END DO
      az=az/(STP-timeinit+1)*10+dcell(3,STP)
        WRITE(22,107) "z","space filling"
  107   FORMAT(A,8X,A13)
      DO divz=1,500
      DO I=timeinit+9,STP,10
       num=0
       NumC=0
       z=dcell(3,I)/500*(divz-0.5d0)
       DO divx=1,100
        x=dcell(1,I)/100*(divx-0.5d0)
        DO divy=1,100
         y=dcell(2,I)/100*(divy-0.5d0)
         NumC=NumC+1
         dist=50
         DO Q=IONS(1,1),IONS(2,1)
          dx=P(1,Q,I)-x
          IF (dx.GE.dcell(1,I)/2d0) THEN
           dx=dx-dcell(1,I)
          ELSE IF (dx.LT.-dcell(1,I)/2d0) THEN
           dx=dx+dcell(1,I)
          END IF
          dy=P(2,Q,I)-y
          IF (dy.GE.dcell(2,I)/2d0) THEN
           dy=dy-dcell(2,I)
          ELSE IF (dy.LT.-dcell(2,I)/2d0) THEN
           dy=dy+dcell(2,I)
          END IF
          dz=P(3,Q,I)-z
          IF (dz.GE.dcell(3,I)/2d0) THEN
           dz=dz-dcell(3,I)
          ELSE IF (dz.LT.-dcell(3,I)/2d0) THEN
           dz=dz+dcell(3,I)
          END IF
          datoms=dx**2d0+dy**2d0+dz**2d0
          IF (datoms.LT.dist) THEN
           dist=datoms
          END IF
         END DO
         IF (dist.GE.1.52d0) THEN
         dist=50
         DO Q=IONS(1,2),IONS(2,2)
          dx=P(1,Q,I)-x
          IF (dx.GE.dcell(1,I)/2d0) THEN
           dx=dx-dcell(1,I)
          ELSE IF (dx.LT.-dcell(1,I)/2d0) THEN
           dx=dx+dcell(1,I)
          END IF
          dy=P(2,Q,I)-y
          IF (dy.GE.dcell(2,I)/2d0) THEN
           dy=dy-dcell(2,I)
          ELSE IF (dy.LT.-dcell(2,I)/2d0) THEN
           dy=dy+dcell(2,I)
          END IF
          dz=P(3,Q,I)-z
          IF (dz.GE.dcell(3,I)/2d0) THEN
           dz=dz-dcell(3,I)
          ELSE IF (dz.LT.-dcell(3,I)/2d0) THEN
           dz=dz+dcell(3,I)
          END IF
          datoms=dx**2d0+dy**2d0+dz**2d0
          IF (datoms.LT.dist) THEN
           dist=datoms
          END IF
         END DO
          IF (datoms.LT.dist) THEN
           dist=datoms
          END IF
         IF (dist.GE.1.20d0) THEN
          num=num+1
         END IF
         END IF
        END DO
       END DO
      END DO
       space=dble(num)/dble(numC)*100d0
       adz=az/500*(divz-0.5d0)
            WRITE(22,108) adz,space
  108   FORMAT(2(F8.5,1X))
      END DO
C
C------close files
        CLOSE(11)
        CLOSE(12)
        CLOSE(22)
        CLOSE(31)
C
        STOP
        END
C


