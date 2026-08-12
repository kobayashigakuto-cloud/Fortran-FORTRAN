C-------declaration of variables--------------
C
        implicit none
        INTEGER Tinit,dr,stepN,Tlast
        INTEGER O,H,ST,Onum1,Onum2
        INTEGER I,J,K
        REAL(8) dx,dy,dz,r,datoms
        REAL(8) Oz1,Oz2,OHdist1,OHdist2
        REAL(8) dvol,z,dOHangle,avevol,avez,OHangle
        REAL(8) OO,Hangle,OHanglelast
        INTEGER num
        CHARACTER(4) ATOM(10)
        INTEGER NTION,NCOMPO,STP,NION(10),IONS(2,10)
        INTEGER IPV(3,5000),preIPV(3,5000)
        REAL(8) CELL(3),P(3,3000,200000),dcell(3,200000)
C-------open files
C
        OPEN(11,FILE='file07.dat',STATUS='old')
        OPEN(12,file='file09p.dat',status='old')
        OPEN(21,file='danglingOH.dat',status='replace')
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
        READ (31,104) Tinit,dr,stepN,Tlast
  104   FORMAT(I8/F8.2/I8/I8)
C          
C-------reading of 'file09p'
C
      DO 210 I=1,Tlast
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
C-------calculate
C
        WRITE(21,107) "Onum","Hnum","averageZ","angle"
  107   FORMAT(A4,2X,A4,4X,A8,8X,A5)

      DO H=IONS(1,2),IONS(2,2)
       num=0
       DO ST=Tinit,Tlast
        OHdist1=10
        OHdist2=20
        DO O=IONS(1,1),IONS(2,1)
         dx=P(1,H,ST)-P(1,O,ST)
         IF (dx.GT.(dcell(1,ST)/2d0)) THEN
          dx=dx-dcell(1,ST)
         END IF
         IF (dx.LT.(-dcell(1,ST)/2d0)) THEN
          dx=dx+dcell(1,ST)
         END IF
         dy=P(2,H,ST)-P(2,O,ST)
         IF (dy.GT.(dcell(2,ST)/2d0)) THEN
          dy=dy-dcell(2,ST)
         END IF
         IF (dy.LT.(-dcell(2,ST)/2d0)) THEN
          dy=dy+dcell(2,ST)
         END IF
         dz=P(3,H,ST)-P(3,O,ST)
         IF (dz.GT.(dcell(3,ST)/2d0)) THEN
          dz=dz-dcell(3,ST)
         END IF
         IF (dz.LT.(-dcell(3,ST)/2d0)) THEN
          dz=dz+dcell(3,ST)
         END IF
         datoms=(dx**2d0+dy**2d0+dz**2d0)**0.5d0
         IF (datoms.LT.OHdist1) THEN
          OHdist2=OHdist1
          Onum2=Onum1
          Oz2=Oz1
          OHdist1=datoms
          Onum1=O
          Oz1=P(3,O,ST)
         ELSE IF (datoms.LT.OHdist2) THEN
          OHdist2=datoms
          Onum2=O
          Oz2=P(3,O,ST)
         ENDIF
        END DO
        IF (OHdist2.LT.2.5) THEN
         num=num+1
        END IF
       END DO
       IF (num/(Tlast-Tinit+1).LT.0.3) THEN
        dvol=0
        z=0
        dOHangle=0
        dx=P(1,H,Tlast)-P(1,Onum1,Tlast)
        IF (dx.GT.(dcell(1,Tlast)/2d0)) THEN
         dx=dx-dcell(1,Tlast)
        END IF
        IF (dx.LT.(-dcell(1,Tlast)/2d0)) THEN
         dx=dx+dcell(1,Tlast)
        END IF
        dy=P(2,H,Tlast)-P(2,Onum1,Tlast)
        IF (dy.GT.(dcell(2,Tlast)/2d0)) THEN
         dy=dy-dcell(2,Tlast)
        END IF
        IF (dy.LT.(-dcell(2,Tlast)/2d0)) THEN
         dy=dy+dcell(2,Tlast)
        END IF
        dz=P(3,H,Tlast)-P(3,Onum1,Tlast)
        IF (dz.GT.(dcell(3,Tlast)/2d0)) THEN
         dz=dz-dcell(3,Tlast)
        END IF
        IF (dz.LT.(-dcell(3,Tlast)/2d0)) THEN
         dz=dz+dcell(3,Tlast)
        END IF
        r=(dx**2+dy**2)**0.5d0
        IF (r.NE.0) THEN
         OHanglelast=atan(dz/r)
        ELSE
         OHanglelast=2*atan(1d0)
        END IF
        DO ST=Tinit,Tlast
         dvol=dvol+dcell(3,ST)-dcell(3,Tlast)
         IF (P(3,Onum1,ST).GE.P(3,Onum1,Tlast)+5) THEN
          z=z+P(3,Onum1,ST)-dcell(3,ST)-P(3,Onum1,Tlast)
         ELSE IF (P(3,Onum1,ST).LT.P(3,Onum1,Tlast)-5) THEN
          z=z+P(3,Onum1,ST)+dcell(3,Tlast)-P(3,Onum1,Tlast)
         ELSE
         z=z+P(3,Onum1,ST)-P(3,Onum1,Tlast)
         END IF
         dx=P(1,H,ST)-P(1,Onum1,ST)
         IF (dx.GT.(dcell(1,ST)/2d0)) THEN
          dx=dx-dcell(1,ST)
         END IF
         IF (dx.LT.(-dcell(1,ST)/2d0)) THEN
          dx=dx+dcell(1,ST)
         END IF
         dy=P(2,H,ST)-P(2,Onum1,ST)
         IF (dy.GT.(dcell(2,ST)/2d0)) THEN
          dy=dy-dcell(2,ST)
         END IF
         IF (dy.LT.(-dcell(2,ST)/2d0)) THEN
          dy=dy+dcell(2,ST)
         END IF
         dz=P(3,H,ST)-P(3,Onum1,ST)
         IF (dz.GT.(dcell(3,ST)/2d0)) THEN
          dz=dz-dcell(3,ST)
         END IF
         IF (dz.LT.(-dcell(3,ST)/2d0)) THEN
          dz=dz+dcell(3,ST)
         END IF
         r=(dx**2+dy**2)**0.5d0
         IF (r.NE.0) THEN
          OHangle=atan(dz/r)
         ELSE
          OHangle=2*atan(1d0)
         END IF
         dOHangle=dOHangle+OHangle-OHanglelast
        END DO
        avevol=dvol/(Tlast-Tinit+1)+dcell(3,Tlast)
        avez=z/(Tlast-Tinit+1)+P(3,Onum1,Tlast)
        OHangle=dOHangle/(Tlast-Tinit+1)+OHanglelast
        OHangle=OHangle*45d0/atan(1d0)
        IF (avez.LT.0) THEN
         avez=avez+avevol
        ELSE IF (avez.GE.avevol) THEN
         avez=avez-avevol
        END IF
        WRITE(21,108) Onum1,H,avez,OHangle,"r"
  108   FORMAT(2(I4,2X),F10.7,2X,F11.7,X,A)
       ELSE
        DO ST=Tinit,Tlast
         dx=P(1,Onum1,ST)-P(1,Onum2,ST)
         IF (dx.GT.(dcell(1,ST)/2d0)) THEN
          dx=dx-dcell(1,ST)
         END IF
         IF (dx.LT.(-dcell(1,ST)/2d0)) THEN
          dx=dx+dcell(1,ST)
         END IF
         dy=P(2,Onum1,ST)-P(2,Onum2,ST)
         IF (dy.GT.(dcell(2,ST)/2d0)) THEN
          dy=dy-dcell(2,ST)
         END IF
         IF (dy.LT.(-dcell(2,ST)/2d0)) THEN
          dy=dy+dcell(2,ST)
         END IF
         dz=P(3,Onum1,ST)-P(3,Onum2,ST)
         IF (dz.GT.(dcell(3,ST)/2d0)) THEN
          dz=dz-dcell(3,ST)
         END IF
         IF (dz.LT.(-dcell(3,ST)/2d0)) THEN
          dz=dz+dcell(3,ST)
         END IF
         OO=(dx**2d0+dy**2d0+dz**2d0)**0.5d0
         Hangle=acos(OHdist1**2d0+OHdist2**2d0-OO**2d0)
     *          /(2d0*OHdist1*OHdist2)+Hangle
        END DO
        Hangle=Hangle/(Tlast-Tinit+1)
        IF (Hangle.LT.22d0/9d0*atan(1d0)) THEN
         dvol=0
         z=0
         dOHangle=0
         dx=P(1,H,Tlast)-P(1,Onum1,Tlast)
         IF (dx.GT.(dcell(1,Tlast)/2d0)) THEN
          dx=dx-dcell(1,Tlast)
         END IF
         IF (dx.LT.(-dcell(1,Tlast)/2d0)) THEN
          dx=dx+dcell(1,Tlast)
         END IF
         dy=P(2,H,Tlast)-P(2,Onum1,Tlast)
         IF (dy.GT.(dcell(2,Tlast)/2d0)) THEN
          dy=dy-dcell(2,Tlast)
         END IF
         IF (dy.LT.(-dcell(2,Tlast)/2d0)) THEN
          dy=dy+dcell(2,Tlast)
         END IF
         dz=P(3,H,Tlast)-P(3,Onum1,Tlast)
         IF (dz.GT.(dcell(3,Tlast)/2d0)) THEN
          dz=dz-dcell(3,Tlast)
         END IF
         IF (dz.LT.(-dcell(3,Tlast)/2d0)) THEN
          dz=dz+dcell(3,Tlast)
         END IF
         r=(dx**2+dy**2)**0.5d0
         IF (r.NE.0) THEN
          OHanglelast=atan(dz/r)
         ELSE
          OHanglelast=2*atan(1d0)
         END IF
         DO ST=Tinit,Tlast
          dvol=dvol+dcell(3,ST)-dcell(3,Tlast)
          IF (P(3,Onum1,ST).GE.P(3,Onum1,Tlast)+5) THEN
           z=z+P(3,Onum1,ST)-dcell(3,ST)-P(3,Onum1,Tlast)
          ELSE IF (P(3,Onum1,ST).LT.P(3,Onum1,Tlast)-5) THEN
           z=z+P(3,Onum1,ST)+dcell(3,Tlast)-P(3,Onum1,Tlast)
          ELSE
          z=z+P(3,Onum1,ST)-P(3,Onum1,Tlast)
          END IF
          dx=P(1,H,ST)-P(1,Onum1,ST)
          IF (dx.GT.(dcell(1,ST)/2d0)) THEN
           dx=dx-dcell(1,ST)
          END IF
          IF (dx.LT.(-dcell(1,ST)/2d0)) THEN
           dx=dx+dcell(1,ST)
          END IF
          dy=P(2,H,ST)-P(2,Onum1,ST)
          IF (dy.GT.(dcell(2,ST)/2d0)) THEN
           dy=dy-dcell(2,ST)
          END IF
          IF (dy.LT.(-dcell(2,ST)/2d0)) THEN
           dy=dy+dcell(2,ST)
          END IF
          dz=P(3,H,ST)-P(3,Onum1,ST)
          IF (dz.GT.(dcell(3,ST)/2d0)) THEN
           dz=dz-dcell(3,ST)
          END IF
          IF (dz.LT.(-dcell(3,ST)/2d0)) THEN
           dz=dz+dcell(3,ST)
          END IF
          r=(dx**2+dy**2)**0.5d0
          IF (r.NE.0) THEN
           OHangle=atan(dz/r)
          ELSE
           OHangle=2*atan(1d0)
          END IF
          dOHangle=dOHangle+OHangle-OHanglelast
         END DO
         avevol=dvol/(Tlast-Tinit+1)+dcell(3,Tlast)
         avez=z/(Tlast-Tinit+1)+P(3,Onum1,Tlast)
         OHangle=dOHangle/(Tlast-Tinit+1)+OHanglelast
         OHangle=OHangle*45d0/atan(1d0)
         IF (avez.LT.0) THEN
          avez=avez+avevol
         ELSE IF (avez.GE.avevol) THEN
          avez=avez-avevol
         END IF
         WRITE(21,108) Onum1,H,avez,OHangle,"a"
        END IF
       END IF
      END DO
C
C------close files
        CLOSE(11)
        CLOSE(12)
        CLOSE(21)
        CLOSE(31)
        END