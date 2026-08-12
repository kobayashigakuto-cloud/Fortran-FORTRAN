C-------declaration of variables--------------
C
        implicit none
        INTEGER time,stepN,STP
        INTEGER O,D
        INTEGER I,J,K,Q,R
        REAL(8) dx,dy,dz,datoms
        REAL(8) dr,dist,dist1,dist2,dist3,dist4,dist5
        INTEGER num
        REAL(8) rate
        CHARACTER(4) ATOM(10)
        INTEGER NTION,NCOMPO,NION(10),IONS(2,10)
        INTEGER IPV(3,5000),preIPV(3,5000)
        REAL(8) CELL(3),P(3,3000,200000),dcell(3,200000)
C
C-------open files
C
        OPEN(11,file='file07.dat',status='old')
        OPEN(12,file='file09p.dat',status='old')
        OPEN(21,file='d5.dat',status='replace')
        OPEN(31,file='parardf.dat',status='old')
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
C-------reading of 'parardf.dat'
C
        READ (31,104) time,dr,stepN,STP
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
C-------rdf
C
        WRITE(21,111) "r (A)","d5_rate"
  111   FORMAT(A5,4X,A7)
       DO dr=2.0,5.0,0.01
       num=0
       DO O=IONS(1,1),IONS(2,1)
        dist=0
        DO I=time,STP
         dist1=50
         dist2=60
         dist3=70
         dist4=80
         dist5=90
         DO D=IONS(1,1),IONS(2,1)
          IF (O.NE.D) THEN
          CALL distance (O,D,I,P,dcell,dx,dy,dz,datoms)
          END IF
          IF (datoms.LT.dist1) THEN
          dist5=dist4
          dist4=dist3
          dist3=dist2
          dist2=dist1
          dist1=datoms
          ELSE IF (datoms.LT.dist2) THEN
          dist5=dist4
          dist4=dist3
          dist3=dist2
          dist2=datoms
          ELSE IF (datoms.LT.dist3) THEN
          dist5=dist4
          dist4=dist3
          dist3=datoms
          ELSE IF (datoms.LT.dist4) THEN
          dist5=dist4
          dist4=datoms
          ELSE IF (datoms.LT.dist5) THEN
          dist5=datoms
          END IF
         END DO
        dist=dist+dist5
        END DO
        IF (dist/(STP-time+1).LT.dr) THEN
        num=num+1
        END IF
       END DO
       rate=dble(num)/dble((IONS(2,1)-IONS(1,1)+1))
        WRITE(21,112) dr,rate
  112   FORMAT(F5.2,4X,F5.2)
       END DO
C
C------close files
        CLOSE(11)
        CLOSE(12)
        CLOSE(21)
        CLOSE(31)
C
        STOP
        END
C
C------subroutine
C
      subroutine distance (Q,R,S,P,dcell,dx,dy,dz,datoms) 
      implicit none
      INTEGER, INTENT(IN)  :: Q,R,S
      REAL(8), INTENT(IN)  :: P(3,3000,125000),dcell(3,125000)
      REAL(8), INTENT(OUT) :: dx,dy,dz,datoms
       dx=P(1,Q,S)-P(1,R,S)
       IF (dx.GT.(dcell(1,S)/2d0)) dx=dx-dcell(1,S)
       IF (dx.LT.(-dcell(1,S)/2d0)) dx=dx+dcell(1,S)
       dy=P(2,Q,S)-P(2,R,S)
       IF (dy.GT.(dcell(2,S)/2d0)) dy=dy-dcell(2,S)
       IF (dy.LT.(-dcell(2,S)/2d0)) dy=dy+dcell(2,S)
       dz=P(3,Q,S)-P(3,R,S)
       IF (dz.GT.(dcell(3,S)/2d0)) dz=dz-dcell(3,S)
       IF (dz.LT.(-dcell(3,S)/2d0)) dz=dz+dcell(3,S)
       datoms=(dx**2d0+dy**2d0+dz**2d0)**0.5d0
      end subroutine
