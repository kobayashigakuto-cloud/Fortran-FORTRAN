C     --- Variables ---
      IMPLICIT NONE
      INTEGER T, CAT, P_IDX, J, J1, J2, CNT
      CHARACTER*30 INFILE
      CHARACTER*30 OUTF(5)
      CHARACTER*20 P_NAMES(5)
C     Buffer: 2 Categories x 10 Values (5 params * 2 val)
      REAL*8 BUF(2, 10)
C
C     --- Define Output Filenames (Per Parameter) ---
      OUTF(1) = '01_OH_dist.dat'
      OUTF(2) = '02_OH_angle.dat'
      OUTF(3) = '03_O_H_dist.dat'
      OUTF(4) = '04_H_angle.dat'
      OUTF(5) = '05_HOH_angle.dat'
C
C     --- Define Parameter Names ---
      P_NAMES(1) = 'OH distance'
      P_NAMES(2) = 'OH angle'
      P_NAMES(3) = 'O..H distance'
      P_NAMES(4) = 'H-bond Angle'
      P_NAMES(5) = 'HOH Angle'
C
C     --- Open Output Files and Write Headers ---
      DO P_IDX = 1, 5
        OPEN(20+P_IDX, FILE=OUTF(P_IDX), STATUS='UNKNOWN')
C
C       Write Legend
C       Col Order: Flat/Dang, Flat/NoDang
        WRITE(20+P_IDX, '(A,A)') 'Summary for: ', P_NAMES(P_IDX)
        WRITE(20+P_IDX, 100) 'Temp(K)',
     &    'Hbond_av', 'Hbond_sd',
     &    'Dang_av',  'Dang_sd'
      END DO
C
  100 FORMAT(A7, 1X, 8(A10, 1X))
C
C     --- Loop over Temperatures (10K to 120K) ---
      DO T = 10, 120, 10
C       Generate filename (e.g. statistics_hb_10K.dat)
C       Adjust format if your filenames are different (e.g. just 10K.dat)
        WRITE(INFILE, 101) T
  101   FORMAT(I0,'K.dat')
C
        OPEN(10, FILE=INFILE, STATUS='OLD')
C
C       Read 4 categories into Buffer
C       Row 1: Flat/Dang, Row 2: Flat/NoDang, etc.
        DO CAT = 1, 2
          READ(10, *) CNT, (BUF(CAT, J), J=1, 10)
        END DO
        CLOSE(10)
C
C       Write Data to each Parameter File
        DO P_IDX = 1, 5
C         Calculate column indices for current parameter
C         P=1 -> J1=1, J2=2 (Avg, Std)
          J1 = (P_IDX - 1) * 2 + 1
          J2 = J1 + 1
C
          WRITE(20+P_IDX, 200) T,
     &      BUF(1, J1), 1.96*BUF(1, J2),
     &      BUF(2, J1), 1.96*BUF(2, J2)
        END DO
      END DO
C
  200 FORMAT(I7, 1X, 8(F10.4, 1X))
C
C     --- Close Output Files ---
      DO P_IDX = 1, 5
        CLOSE(20+P_IDX)
      END DO
C
      END