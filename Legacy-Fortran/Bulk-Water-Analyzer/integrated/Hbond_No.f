C-------declaration of variables--------------
C
        implicit none
      INTEGER I,J,K,ST,Tinit,Tlast,STP
      INTEGER H,O,Onum1,Onum2,Fst,Snd
      REAL(8) Fstdist,Snddist,Hangle,PI
      CHARACTER(4) ATOM(10)
      INTEGER NTION,NCOMPO,NION(10),IONS(2,10)
      INTEGER IPV(3,3000)
      REAL(8) CELL(3),P(3,3000,125000),dcell(3,125000)
      REAL(8) dr
      INTEGER stepN
      INTEGER MAX_O
      PARAMETER (MAX_O=3000)
      INTEGER HB_COUNT(MAX_O)
      REAL(8) HB_HIST(0:6)
      PI = ACOS(-1.0d0)
C
C-------open files-------
C
      OPEN(11,FILE='file07.dat',STATUS='old')
      OPEN(12,FILE='file09p.dat',STATUS='old')
      OPEN(31,FILE='parardf.dat',STATUS='old')
      OPEN(41,FILE='hb_distribution.dat',STATUS='replace')
C
C-------reading of 'file07'-------
C
      READ (11,101) NTION,NCOMPO
      READ (11,102) (ATOM(I),I=1,NCOMPO)
      READ (11,103) (NION(I),I=1,NCOMPO)
      READ (11,103) (IONS(1,I),I=1,NCOMPO)
      READ (11,103) (IONS(2,I),I=1,NCOMPO)
  101 FORMAT (/I7,I3)
  102 FORMAT (10(2X,A4))
  103 FORMAT (10I6)

C-------reading of 'parardf.dat'-------
      READ (31,104) Tinit,dr,stepN,Tlast
  104 FORMAT(I8/F8.2/I8/I8)

C-------reading of 'file09p'-------
      DO I=1,Tlast
        READ (12,105) (CELL(K),K=1,3)
        READ (12,106) ((IPV(K,J),K=1,3),J=1,NTION)
  105   FORMAT (10X,2(F7.3,21X),F7.3)
  106   FORMAT (18I4)
        dcell(1,I)=CELL(1)
        dcell(2,I)=CELL(2)
        dcell(3,I)=CELL(3)
        DO J=1,NTION
          DO K=1,3
            P(K,J,I)=CELL(K)*DBLE(IPV(K,J))/9d3
          END DO
        END DO
      END DO
C
C-------main loop-------
C
      DO I=0,6
        HB_HIST(I)=0.0d0
      END DO

      DO ST=Tinit,Tlast
        DO J=1,NTION
          HB_COUNT(J)=0
        END DO

        DO H=IONS(1,2),IONS(2,2) 
          CALL base (1,H,ST,P,dcell,IONS,Onum1,Onum2,Fstdist,Snddist)

          CALL triangle (ST,ST,H,Onum1,Onum2,P,dcell,Hangle)
          Hangle=PI-Hangle 

          IF (Snddist.LT.2.5d0.AND.Hangle.LT.(7d0/36d0*PI)) THEN
            HB_COUNT(Onum1)=HB_COUNT(Onum1)+1
            HB_COUNT(Onum2)=HB_COUNT(Onum2)+1
          END IF
        END DO

        DO O=361,720
          J=HB_COUNT(O)
          IF (J.GT.6) J=6
          HB_HIST(J)=HB_HIST(J)+1.0d0
        END DO
      END DO

      DO I=0,6
        WRITE(41,'(I,F20.8)') I,HB_HIST(I)/((Tlast-Tinit+1)*360.0d0)
      END DO
C
C-------close files-------
      CLOSE(11)
      CLOSE(12)
      CLOSE(31)
      CLOSE(41)
      
      STOP
      END PROGRAM

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

      subroutine arccos (OH1,OH2,HH,dHOH) 
        implicit none
        REAL(8), INTENT(IN) :: OH1,OH2,HH
        REAL(8), INTENT(INOUT) :: dHOH
        REAL(8) :: acos_val
        IF (OH1*OH2.GT.1.0d-6) THEN
          acos_val=(OH1**2d0+OH2**2d0-HH**2d0)/(2d0*OH1*OH2)
          IF (acos_val.GT. 1.0d0) acos_val= 1.0d0
          IF (acos_val.LT.-1.0d0) acos_val=-1.0d0
          dHOH=acos(acos_val)+dHOH
        END IF
      end subroutine

      subroutine base (a,Onum,T,P,dcell,IONS,Fst,Snd,Fstdist,Snddist) 
      implicit none
      INTEGER, INTENT(IN)  :: T,Onum,a
      REAL(8), INTENT(IN)  :: P(3,3000,125000),dcell(3,125000)
      INTEGER, INTENT(OUT) :: Fst,Snd
      INTEGER :: Q,IONS(2,10)
      REAL(8) :: Fstdist,Snddist,dx,dy,dz,datoms
        Fstdist=100.0d0
        Snddist=200.0d0
        DO Q=IONS(1,a),IONS(2,a)
        CALL distance (Q,Onum,T,P,dcell,dx,dy,dz,datoms)
          IF (datoms.LT.Fstdist) THEN
           Snddist=Fstdist
           Snd=Fst
           Fstdist=datoms
           Fst=Q
          ELSE IF (datoms.LT.Snddist) THEN
           Snddist=datoms
           Snd=Q
          ENDIF
         END DO
      end subroutine

      subroutine triangle (Tinit,Tlast,A,B,C,P,dcell,BAC) 
        implicit none
        INTEGER, INTENT(IN)  :: Tinit,Tlast,A,B,C
        REAL(8), INTENT(IN)  :: P(3,3000,125000),dcell(3,125000)
        REAL(8), INTENT(OUT) :: BAC
        INTEGER :: ST
        REAL(8) :: d,AB,AC,BC,dx,dy,dz
       d=0.0d0
       DO ST=Tinit,Tlast
       CALL distance (A,B,ST,P,dcell,dx,dy,dz,AB)
       CALL distance (A,C,ST,P,dcell,dx,dy,dz,AC)
       CALL distance (B,C,ST,P,dcell,dx,dy,dz,BC)
       CALL arccos (AB,AC,BC,d)
         END DO
         BAC=d/DBLE(Tlast-Tinit+1)
      end subroutine