C-------declaration of variables--------------
C
        implicit none
        INTEGER time,stepN,STP
        INTEGER O,D
        INTEGER I,J,K
        REAL(8) dx,dy,dz,datoms
        REAL(8) dr,dist1,dist2,dist3,dist4,dist5
        REAL(8) dist_1
        REAL(8) scale_val
        INTEGER num1
        INTEGER NumA1,NumB1
        REAL(8) rate
        
        REAL(8) val_d5(5000)
        INTEGER total_steps

        CHARACTER(4) ATOM(10)
        INTEGER NTION,NCOMPO,NION(10),IONS(2,10)
        INTEGER IPV(3,5000)
        REAL(8) CELL(3),P(3,3000,200000),dcell(3,200000)
C
C-------open files
C
        OPEN(11,file='file07.dat',status='old')
        OPEN(12,file='file09p.dat',status='old')
        OPEN(21,file='d5_ratio_bulk.dat',status='replace')
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
      DO I=1,STP
        READ (12,105) (CELL(K),K=1,3)
        READ (12,106) ((IPV(K,J),K=1,3),J=1,NTION)
  105   FORMAT (10X,F7.3,21X,F7.3,21X,F7.3)
  106   FORMAT (18I4)
            dcell(1,I)=CELL(1)
            dcell(2,I)=CELL(2)
            dcell(3,I)=CELL(3)
           DO J=1,NTION
           DO K=1,3
               P(K,J,I)=CELL(K)*dble(IPV(K,J))/9d3
           END DO
           END DO
      END DO
C
C-------distance
C
       total_steps = STP - time + 1

       DO O=IONS(1,1),IONS(2,1)
        num1=0; dist_1=0.0d0
        
        DO I=time,STP

         dist1=50.0d0; dist2=60.0d0; dist3=70.0d0
         dist4=80.0d0; dist5=90.0d0
         
         DO D=IONS(1,1),IONS(2,1)
          IF (O.NE.D) THEN
           CALL distance (O,D,I,P,dcell,dx,dy,dz,datoms)
          ELSE
           datoms=1000.0d0
          END IF
          
          IF (datoms.LT.dist1) THEN
             dist5=dist4; dist4=dist3; dist3=dist2; dist2=dist1; dist1=datoms
          ELSE IF (datoms.LT.dist2) THEN
             dist5=dist4; dist4=dist3; dist3=dist2; dist2=datoms
          ELSE IF (datoms.LT.dist3) THEN
             dist5=dist4; dist4=dist3; dist3=datoms
          ELSE IF (datoms.LT.dist4) THEN
             dist5=dist4; dist4=datoms
          ELSE IF (datoms.LT.dist5) THEN
             dist5=datoms
          END IF
         END DO
         
        scale_val = (dist1 + dist2 + dist3) / 3.0d0
        
        num1=num1+1
        dist_1=dist_1+(dist5/scale_val)
        END DO
        
        val_d5(O)=dist_1/dble(num1)
       END DO

       WRITE(21,111) "Ratio","bulk"
  111  FORMAT(A5,4X,A8)

       DO dr=1.00,2.50,0.01
       NumA1=0; NumB1=0
       
       DO O=IONS(1,1),IONS(2,1)
         NumB1=NumB1+1
         IF (val_d5(O).LT.dr) NumA1=NumA1+1
       END DO
       
       IF (NumB1.GT.0) THEN
        rate=dble(NumA1)/dble(NumB1)
        WRITE(21,112) dr,rate
       END IF
       END DO
  112  FORMAT(F5.2,4X,F10.5)

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
      REAL(8), INTENT(IN)  :: P(3,3000,200000),dcell(3,200000)
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