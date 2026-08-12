C-------declaration of variables-------
C
        implicit none
        INTEGER I,J,K,Tinit,dr,stepN,Tlast
        INTEGER H,O,ST,Onum1,Onum2
        INTEGER num,Fst,Snd
        REAL(8) dx,dy,dz,datoms,r,acos_val,vol,z,Fstdist,Snddist
        REAL(8) OHdist,O_Hdist
        REAL(8) Hangle,HOH
        REAL(8) OHangle,dOHangle
        REAL(8) unitOH(3)
        CHARACTER(4) ATOM(10)
        INTEGER NTION,NCOMPO,STP,NION(10),IONS(2,10)
        INTEGER IPV(3,3000),preIPV(3,3000)
        REAL(8) CELL(3),P(3,3000,125000),dcell(3,125000)
C
C-------Added variables for statistics-------
C       CAT: 1=H-bond, 2=Dangling
C       VALS: Stores [z, OHdist, OHangle, O_Hdist, Hangle, HOH] (6 items)
        INTEGER S_CNT(2), CAT, M
        REAL(8) S_SUM(2,6), S_SQ(2,6), VALS(6)
        REAL(8) AVG_VAL, VAR_VAL, STD_VAL
C
C-------open files
C
        OPEN(11,file='file07.dat',status='old')
        OPEN(12,file='file09p.dat',status='old')
        OPEN(21,file='Parameter_Hbond.dat',status='replace')
        OPEN(22,file='Parameter_dangling.dat',status='replace')
C       New statistics file
        OPEN(25,file='statistics_hb_2cat.dat',status='replace')
        OPEN(31,file='parardf.dat',status='old')
C
C-------Initialize Statistics Arrays
C
        DO I=1,2
          S_CNT(I) = 0
          DO J=1,6
            S_SUM(I,J) = 0.0d0
            S_SQ(I,J)  = 0.0d0
          END DO
        END DO
C
C-------reading of 'file07','parardf.dat'
C
        READ (11,101) NTION,NCOMPO
        READ (11,102) (ATOM(I),I=1,NCOMPO)
        READ (11,103) (NION(I),I=1,NCOMPO)
        READ (11,103) (IONS(1,I),I=1,NCOMPO)
        READ (11,103) (IONS(2,I),I=1,NCOMPO)
  101   FORMAT (/I7,I3)
  102   FORMAT (10(2X,A4))
  103   FORMAT (10I6)
        READ (31,104) Tinit,dr,stepN,Tlast
  104   FORMAT(I8/F8.2/I8/I8)
C
C-------reading of 'file09p'
C
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
          P(K,J,I)=CELL(K)*dble(IPV(K,J))/9d3
          END DO
        END DO
      END DO
C
C-------index
C
        WRITE(21,107) "Hnum","averageZ","HOdist","HOangle",
     *  "OHdist","OHangle","HOHangle"
        WRITE(22,107) "Hnum","averageZ","HOdist","HOangle",
     *  "OHdist","OHangle","HOHangle"
  107   FORMAT(A4,X,6A12)
C
C--------calculate----------
C
      vol=0
      DO ST=Tinit,Tlast
       vol=vol+dcell(3,ST)
      END DO
      vol=vol/(Tlast-Tinit+1)

      DO H=IONS(1,2),IONS(2,2)
       num=0
       OHdist=0d0; O_Hdist=0d0
       DO ST=Tinit,Tlast
       CALL base (1,H,ST,P,dcell,IONS,Onum1,Onum2,Fstdist,Snddist)
        OHdist =Fstdist+OHdist
        O_Hdist=Snddist+OHdist
        IF (Snddist.LT.2.5) num=num+1
       END DO
       OHdist =OHdist /(Tlast-Tinit+1)
       O_Hdist=O_Hdist/(Tlast-Tinit+1)

       CALL triangle (Tinit,Tlast,H,Onum1,Onum2,P,dcell,Hangle)
        Hangle=acos(-1d0)-Hangle
        
C      Judge H-bond
       IF ((REAL(num,8)/REAL(Tlast-Tinit+1,8).GE.0.9d0).AND.
     *     (Hangle.LT.7d0/36d0*acos(-1d0))) THEN
C        --- Case: H-bond (Bonded) ---
         z=0
         DO ST=Tinit,Tlast
          z=z+P(3,H,ST)
         END DO
         z=z/(Tlast-Tinit+1)
         IF (z.LT.0) THEN
          z=z+vol
         ELSE IF (z.GE.vol) THEN
          z=z-vol
         END IF

         dOHangle=0
         DO ST=Tinit,Tlast
       CALL distance (H,Onum1,ST,P,dcell,dx,dy,dz,datoms)
          unitOH(1)=dx/datoms
          unitOH(2)=dy/datoms
          unitOH(3)=dz/datoms
          r=(unitOH(1)**2+unitOH(2)**2)**0.5d0
          IF (r.NE.0) THEN
           OHangle=atan(unitOH(3)/r)
          ELSE
           OHangle=2*atan(1d0)
          END IF
          dOHangle=dOHangle+OHangle
         END DO
         OHangle=dOHangle/(Tlast-Tinit+1)
         OHangle=OHangle*180d0/acos(-1d0)
         Hangle=acos(-1d0)-Hangle
       CALL base (2,Onum1,Tlast,P,dcell,IONS,Fst,Snd,Fstdist,Snddist)
       CALL triangle (Tinit,Tlast,Onum1,Fst,Snd,P,dcell,HOH)
         HOH=HOH*180d0/acos(-1d0)
         
         WRITE(21,108) H,z,OHdist,OHangle,O_Hdist,Hangle,HOH
         CAT = 1  ! H-bond Category

       ELSE
C        --- Case: No H-bond (Dangling) ---
         z=0
         DO ST=Tinit,Tlast
          z=z+P(3,H,ST)
         END DO
         z=z/(Tlast-Tinit+1)
         IF (z.LT.0) THEN
          z=z+vol
         ELSE IF (z.GE.vol) THEN
          z=z-vol
         END IF

         dOHangle=0
         DO ST=Tinit,Tlast
       CALL distance (H,Onum1,ST,P,dcell,dx,dy,dz,datoms)
          unitOH(1)=dx/datoms
          unitOH(2)=dy/datoms
          unitOH(3)=dz/datoms
          r=(unitOH(1)**2+unitOH(2)**2)**0.5d0
          IF (r.NE.0) THEN
           OHangle=atan(unitOH(3)/r)
          ELSE
           OHangle=2*atan(1d0)
          END IF
          dOHangle=dOHangle+OHangle
         END DO
         OHangle=dOHangle/(Tlast-Tinit+1)
         OHangle=OHangle*180d0/acos(-1d0)
         Hangle=acos(-1d0)-Hangle
       CALL base (2,Onum1,Tlast,P,dcell,IONS,Fst,Snd,Fstdist,Snddist)
       CALL triangle (Tinit,Tlast,Onum1,Fst,Snd,P,dcell,HOH)
         HOH=HOH*180d0/acos(-1d0)
         
         WRITE(22,108) H,z,OHdist,OHangle,O_Hdist,Hangle,HOH
         CAT = 2  ! Dangling Category
       END IF
       
  108   FORMAT(I4,X,10F12.7)

C      --- Accumulate Statistics ---
C      Store values in temp array [z, OHdist, OHangle, O_Hdist, Hangle, HOH]
       VALS(1) = z
       VALS(2) = OHdist
       VALS(3) = OHangle
       VALS(4) = O_Hdist
       VALS(5) = Hangle
       VALS(6) = HOH
       
       S_CNT(CAT) = S_CNT(CAT) + 1
       DO M=1,6
         S_SUM(CAT,M) = S_SUM(CAT,M) + VALS(M)
         S_SQ(CAT,M)  = S_SQ(CAT,M) + VALS(M)**2d0
       END DO

      END DO

C
      DO cat = 1, 2
        WRITE(25, '(I8)', advance='no') s_cnt(cat)
        DO m = 2, 6
            IF (s_cnt(cat) .GT. 1) THEN
             avg_val = s_sum(cat, m) / REAL(s_cnt(cat), 8)
             ! Variance = (SumSq - Sum^2/N) / (N-1)
             var_val = (s_sq(cat,m) - (s_sum(cat,m)**2d0)/
     *                  REAL(s_cnt(cat),8)) / REAL(s_cnt(cat)-1,8)
             IF (var_val .LT. 0d0) var_val = 0d0
             std_val = SQRT(var_val)
            ELSE IF (s_cnt(cat) .EQ. 1) THEN
             avg_val = s_sum(cat, m)
             std_val = 0.0d0
            ELSE
             avg_val = 0.0d0
             std_val = 0.0d0
            END IF
            WRITE(25, '(2F12.4)', advance='no') avg_val, std_val
        END DO
        WRITE(25, *) ! New line
      END DO

C
C-------close files
C
        CLOSE(11)
        CLOSE(12)
        CLOSE(21)
        CLOSE(22)
        CLOSE(25)
        CLOSE(31)
        END
C
C-------subroutine
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
        Fstdist=100
        Snddist=200
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
       d=0
       DO ST=Tinit,Tlast
       CALL distance (A,B,ST,P,dcell,dx,dy,dz,AB)
       CALL distance (A,C,ST,P,dcell,dx,dy,dz,AC)
       CALL distance (B,C,ST,P,dcell,dx,dy,dz,BC)
       CALL arccos (AB,AC,BC,d)
         END DO
         BAC=d/(Tlast-Tinit+1)
      end subroutine