C-------declaration of variables--------------
C
        implicit none
        INTEGER TNST,NST,stepN,timeinit,DST,num
        INTEGER O,H
        INTEGER I,J,K,KM,Q
        REAL(8) dx,dy,dz,r,theta,rate
        REAL(8) angle(200)
        REAL(8) dvol,avevol,datoms,dr
        INTEGER Nstart,dstep
        INTEGER nO1H,nO2H
        REAL(8) O1H,O2H,dO1H,dO2H
        REAL(8) ax,ay,az,dOHx,dOHy,dOHz
        REAL(8) O1Hx,O1Hy,O1Hz,O2Hx,O2Hy,O2Hz
        CHARACTER(4) ATOM(10)
        INTEGER NTION,NCOMPO,STP,NION(10),IONS(2,10)
        INTEGER IPV(3,5000),preIPV(3,5000)
        REAL(8) CELL(3),P(3,3000,200000),dcell(3,200000)
C-------open files
C
        OPEN(11,FILE='file07.dat',STATUS='old')
        OPEN(12,file='file09p.dat',status='old')
        OPEN(22,file='modified.dat',status='replace')
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
	timeinit=STP-(STP-timeinit+1)*10+1
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
C-------Calculate O-H distance
C
	DO K = 10,180,10
	 angle(K)=0d0
	END DO
	num = 0
       DO Q = IONS(1,1),IONS(2,1)
        DO I = timeinit,STP
	 num=num+1
         O1H = 100d0
         O2H = 100d0
         nO1H = 1000d0
         nO2H = 1000d0
         DO J = IONS(1,2),IONS(2,2)
            dx=P(1,J,I)-P(1,Q,I)
            IF (dx.GT.(dcell(1,I)/2d0)) THEN
               dx=dx-dcell(1,I)
            END IF
            IF (dx.LT.(-dcell(1,I)/2d0)) THEN
               dx=dx+dcell(1,I)
            END IF
               dy=P(2,J,I)-P(2,Q,I)
            IF (dy.GT.(dcell(2,I)/2d0)) THEN
               dy=dy-dcell(2,I)
            END IF
            IF (dy.LT.(-dcell(2,I)/2d0)) THEN
               dy=dy+dcell(2,I)
            END IF
               dz=P(3,J,I)-P(3,Q,I)
            IF (dz.GT.(dcell(3,I)/2d0)) THEN
               dz=dz-dcell(3,I)
            END IF
            IF (dz.LT.(-dcell(3,I)/2d0)) THEN
               dz=dz+dcell(3,I)
            END IF
            datoms=(dx**2d0+dy**2d0+dz**2d0)**0.5d0
            IF (O1H.GE.datoms) THEN
             O2H = O1H
             nO2H = nO1H
             O1H = datoms
             nO1H = J
            ELSE IF (O2H.GE.datoms) THEN
             O2H = datoms
             nO2H = J
            END IF
         END DO
	 O1Hx=P(1,nO1H,I)
	 O1Hy=P(2,nO1H,I)
	 O1Hz=P(3,nO1H,I)
	 O2Hx=P(1,nO2H,I)
	 O2Hy=P(2,nO2H,I)
	 O2Hz=P(3,nO2H,I)	 
         ax = (O1Hx+O2Hx)/2
         ay = (O1Hy+O2Hy)/2
         az = (O1Hz+O2Hz)/2
	 IF (O1Hx-O2Hx.GT.(dcell(1,I)/2d0)) THEN
	  IF (ax.GE.dcell(1,I)/2d0) THEN
	   ax = ax - dcell(1,I)/2d0
	  ELSE
	   ax = ax + dcell(1,I)/2d0
	  END IF
	 ELSE IF (O2Hx-O1Hx.GT.(dcell(1,I)/2d0)) THEN
	  IF (ax.GE.dcell(1,I)/2d0) THEN
	   ax = ax - dcell(1,I)/2d0
	  ELSE
	   ax = ax + dcell(1,I)/2d0
	  END IF
	 END IF
	 IF (O1Hy-O2Hy.GT.(dcell(2,I)/2d0)) THEN
	  IF (ay.GE.dcell(2,I)/2d0) THEN
	   ay = ay - dcell(2,I)/2d0
	  ELSE
	   ay = ay + dcell(2,I)/2d0
	  END IF
	 ELSE IF (O2Hy-O1Hy.GT.(dcell(2,I)/2d0)) THEN
	  IF (ay.GE.dcell(2,I)/2d0) THEN
	   ay = ay - dcell(2,I)/2d0
	  ELSE
	   ay = ay + dcell(2,I)/2d0
	  END IF
	 END IF
	 IF (O1Hz-O2Hz.GT.(dcell(3,I)/2d0)) THEN
	  IF (az.GE.dcell(3,I)/2d0) THEN
	   az = az - dcell(3,I)/2d0
	  ELSE
	   az = az + dcell(3,I)/2d0
	  END IF
	 ELSE IF (O2Hz-O1Hz.GT.(dcell(3,I)/2d0)) THEN
	  IF (az.GE.dcell(3,I)/2d0) THEN
	   az = az - dcell(3,I)/2d0
	  ELSE
	   az = az + dcell(3,I)/2d0
	  END IF
	 END IF
	 dOHx=P(1,Q,I)-ax
	 dOHy=P(2,Q,I)-ay
	 dOHz=P(3,Q,I)-az
	 IF (dOHx.GT.(dcell(1,I)/2d0)) THEN
	  dOHx = dOHx - dcell(1,I)
	 ELSE IF (dOHx.LT.(-dcell(1,I)/2d0)) THEN
	  dOHx = dOHx + dcell(1,I)
	 END IF
	 IF (dOHy.GT.(dcell(2,I)/2d0)) THEN
	  dOHy = dOHy - dcell(2,I)
	 ELSE IF (dOHy.LT.(-dcell(2,I)/2d0)) THEN
	  dOHy = dOHy + dcell(2,I)
	 END IF
	 IF (dOHz.GT.(dcell(3,I)/2d0)) THEN
	  dOHz = dOHz - dcell(3,I)
	 ELSE IF (dOHz.LT.(-dcell(3,I)/2d0)) THEN
	  dOHz = dOHz + dcell(3,I)
	 END IF
        r=(dOHx**2+dOHy**2)**0.5d0
	IF (r.NE.0) THEN
	 theta=atan(dOHz/r)
	theta=2*atan(1d0)-theta
	rate=1/sin(theta)
	IF (rate.GE.1000) THEN
	 rate=1000
	END IF
	ELSE
	 IF (dOHz.GE.0) THEN
	  theta=0
	 ELSE
	  theta=4d0*atan(1d0)
	 END IF
	 rate=1000
	END IF
	theta=theta*45d0/atan(1d0)
	DO K = 10,180,10
	 IF ((theta.GE.K-10).AND.(theta.LT.K)) THEN
	  angle(K)=angle(K)+rate
	 END IF
	END DO
        END DO
       END DO
	DO K = 10,180,10
            WRITE(22,109)angle(K)/num
  109   FORMAT(F15.6)
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