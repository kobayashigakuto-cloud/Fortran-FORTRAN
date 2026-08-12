module mod_constants
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: PI = 3.14159265358979323846_dp
end module mod_constants

module mod_geometry
    use mod_constants
    implicit none
contains
    ! 周期境界条件を考慮した距離計算
    subroutine get_dist_pbc(p1, p2, cell, dist, diff_vec)
        real(dp), intent(in) :: p1(3), p2(3), cell(3)
        real(dp), intent(out) :: dist
        real(dp), intent(out), optional :: diff_vec(3)
        real(dp) :: d(3)
        d = p1 - p2
        d = d - nint(d / cell) * cell
        dist = sqrt(sum(d**2))
        if (present(diff_vec)) diff_vec = d
    end subroutine

    ! 3点間の角度 (BAC) を計算
    subroutine get_angle_3pt(pa, pb, pc, cell, angle_rad)
        real(dp), intent(in) :: pa(3), pb(3), pc(3), cell(3)
        real(dp), intent(out) :: angle_rad
        real(dp) :: v_ab(3), v_ac(3), d_ab, d_ac, cos_t
        call get_dist_pbc(pb, pa, cell, d_ab, v_ab)
        call get_dist_pbc(pc, pa, cell, d_ac, v_ac)
        if (d_ab * d_ac < 1.0e-10_dp) then
            angle_rad = 0.0_dp
        else
            cos_t = sum(v_ab * v_ac) / (d_ab * d_ac)
            angle_rad = acos(max(-1.0_dp, min(1.0_dp, cos_t)))
        end if
    end subroutine
end module mod_geometry

program water_system_analyzer
    use mod_constants
    use mod_geometry
    implicit none

    ! --- 変数宣言 ---
    integer :: n_total, n_comp, t_init, t_last, step_n, t, i, j, k, s, step_count
    integer :: n_o, n_h, i_h, i_o, j_o, fst, snd, num, Onum1, Onum2, cat
    real(dp) :: dr_val, d1, d2, h_ang, datoms, OHdist1, OHdist2, OO, r_threshold
    real(dp) :: cur_z, cur_r, cur_ang, h_angle_avg, hoh_angle, avg_val, std_val, var_val
    
    integer, allocatable :: n_ion_type(:), ions(:,:), ipv(:,:)
    character(4), allocatable :: atom_names(:)
    real(dp), allocatable :: coords(:,:,:), cell_data(:,:)
    
    ! 解析用配列
    real(dp), allocatable :: avg_d5_dist(:)
    real(dp) :: hb_coord_distri(0:6)
    integer, allocatable :: current_hb_count(:)
    
    ! 統計用配列 [Category: 1=HB, 2=Dangling] [Param: 1:z, 2:OHdist, 3:OHangle, 4:O_Hdist, 5:Hangle, 6:HOH]
    integer :: s_cnt(2)
    real(dp) :: s_sum(2, 6), s_sq(2, 6), vals(6)

    ! --- データ入力 ---
    open(11, file='file07.dat', status='old')
    read(11, '(/I7, I3)') n_total, n_comp
    allocate(atom_names(n_comp), n_ion_type(n_comp), ions(2, n_comp))
    read(11, '(10(2X,A4))') atom_names
    read(11, '(10I6)') n_ion_type
    read(11, '(10I6)') ions(1, :); read(11, '(10I6)') ions(2, :)
    close(11)

    open(31, file='parardf.dat', status='old')
    read(31, *) t_init, dr_val, step_n, t_last
    close(31)

    step_count = t_last - t_init + 1
    n_o = ions(2,1) - ions(1,1) + 1
    n_h = ions(2,2) - ions(1,2) + 1

    allocate(coords(3, n_total, t_last), cell_data(3, t_last), ipv(3, n_total))
    
    print *, "Reading Trajectory..."
    open(12, file='file09p.dat', status='old')
    do t = 1, t_last
        read(12, '(10X, F7.3, 21X, F7.3, 21X, F7.3)') cell_data(:, t)
        read(12, '(18I4)') ((ipv(k, j), k=1,3), j=1, n_total)
        do j = 1, n_total
            coords(:, j, t) = cell_data(:, t) * real(ipv(:, j), dp) / 9000.0_dp
        end do
    end do
    close(12)

    ! --- 初期化 ---
    allocate(avg_d5_dist(ions(1,1):ions(2,1))); avg_d5_dist = 0.0_dp
    allocate(current_hb_count(ions(1,1):ions(2,1)))
    hb_coord_distri = 0.0_dp
    s_cnt = 0; s_sum = 0.0_dp; s_sq = 0.0_dp

    open(21, file='Parameter_Hbond.dat', status='replace')
    open(22, file='Parameter_dangling.dat', status='replace')
    open(23, file='d5.dat', status='replace')
    open(24, file='hb_distribution.dat', status='replace')
    open(25, file='statistics_hb_2cat.dat', status='replace')

    write(21, '(A4, 6A12)') "Hnum","averageZ","HOdist","HOangle","OHdist","OHangle","HOHangle"
    write(22, '(A4, 6A12)') "Hnum","averageZ","HOdist","HOangle","OHdist","OHangle","HOHangle"

    ! --- メイン解析ループ1: H原子ベース (HB/Dangling判定) ---
    print *, "Analyzing H-atoms (HB vs Dangling)..."
    do i_h = ions(1, 2), ions(2, 2)
        num = 0
        h_angle_avg = 0.0_dp
        OHdist1 = 0.0_dp; OHdist2 = 0.0_dp
        
        ! ステップループ：最近接酸素の特定と幾何条件の集計
        do t = t_init, t_last
            d1 = 100.0_dp; d2 = 200.0_dp
            do i_o = ions(1, 1), ions(2, 1)
                call get_dist_pbc(coords(:,i_h,t), coords(:,i_o,t), cell_data(:,t), datoms)
                if (datoms < d1) then
                    d2 = d1; snd = fst; d1 = datoms; fst = i_o
                else if (datoms < d2) then
                    d2 = datoms; snd = i_o
                end if
            end do
            if (t == t_last) then; Onum1 = fst; Onum2 = snd; end if
            
            OHdist1 = OHdist1 + d1
            OHdist2 = OHdist2 + d2
            if (d2 < 2.5_dp) num = num + 1
            
            ! ∠O-H...O 角度の積算
            call get_dist_pbc(coords(:,fst,t), coords(:,snd,t), cell_data(:,t), OO)
            ! 余弦定理 (正しい括弧位置)
            h_ang = acos(max(-1.0_dp, min(1.0_dp, (d1**2 + d2**2 - OO**2)/(2.0_dp*d1*d2))))
            h_angle_avg = h_angle_avg + (PI - h_ang)
        end do
        
        h_angle_avg = h_angle_avg / real(step_count, dp)
        OHdist1 = OHdist1 / real(step_count, dp)
        OHdist2 = OHdist2 / real(step_count, dp)

        ! 判定 (num割合 >= 0.9 かつ 角度 < 35度)
        if ((real(num,dp)/real(step_count,dp) >= 0.9_dp) .and. (h_angle_avg < (35.0_dp*PI/180.0_dp))) then
            cat = 1 ! HB
        else
            cat = 2 ! Dangling
        end if

        ! Z座標平均とOH角度平均の計算
        cur_z = 0.0_dp; cur_ang = 0.0_dp
        do t = t_init, t_last
            cur_z = cur_z + coords(3, i_h, t)
            block
                real(dp) :: dv(3), dr
                call get_dist_pbc(coords(:,i_h,t), coords(:,Onum1,t), cell_data(:,t), datoms, dv)
                dr = sqrt(dv(1)**2 + dv(2)**2)
                if (dr > 1.0e-10_dp) then
                    cur_ang = cur_ang + atan(dv(3)/dr)
                else
                    cur_ang = cur_ang + sign(PI/2.0_dp, dv(3))
                end if
            end block
        end do
        cur_z = cur_z / real(step_count, dp)
        cur_ang = (cur_ang / real(step_count, dp)) * 180.0_dp / PI
        
        ! 水分子自体のHOH角度 (Onum1を酸素とする)
        hoh_angle = 0.0_dp
        block
            integer :: o_fst, o_snd
            ! 酸素Onum1に結合している2つの水素を探す (簡易的に最近接2つ)
            d1 = 100.0_dp; d2 = 200.0_dp
            do j = ions(1,2), ions(2,2)
                call get_dist_pbc(coords(:,Onum1,t_last), coords(:,j,t_last), cell_data(:,t_last), datoms)
                if (datoms < d1) then
                    d2 = d1; o_snd = o_fst; d1 = datoms; o_fst = j
                else if (datoms < d2) then
                    d2 = datoms; o_snd = j
                end if
            end do
            call get_angle_3pt(coords(:,Onum1,t_last), coords(:,o_fst,t_last), coords(:,o_snd,t_last), &
                               cell_data(:,t_last), hoh_angle)
            hoh_angle = hoh_angle * 180.0_dp / PI
        end block

        ! 結果の書き出し
        if (cat == 1) then
            write(21, '(I4, X, 10F12.7)') i_h, cur_z, OHdist1, cur_ang, OHdist2, h_angle_avg*180.0_dp/PI, hoh_angle
        else
            write(22, '(I4, X, 10F12.7)') i_h, cur_z, OHdist1, cur_ang, OHdist2, h_angle_avg*180.0_dp/PI, hoh_angle
        end if

        ! 統計積算
        vals = [cur_z, OHdist1, cur_ang, OHdist2, h_angle_avg*180.0_dp/PI, hoh_angle]
        s_cnt(cat) = s_cnt(cat) + 1
        s_sum(cat, :) = s_sum(cat, :) + vals
        s_sq(cat, :) = s_sq(cat, :) + vals**2
    end do

    ! --- 統計出力 ---
    do cat = 1, 2
        write(25, '(I8)', advance='no') s_cnt(cat)
        do j = 1, 6
            if (s_cnt(cat) > 1) then
                avg_val = s_sum(cat, j) / real(s_cnt(cat), dp)
                var_val = (s_sq(cat, j) - (s_sum(cat, j)**2)/real(s_cnt(cat), dp)) / real(s_cnt(cat)-1, dp)
                std_val = sqrt(max(0.0_dp, var_val))
            else
                avg_val = merge(s_sum(cat, j), 0.0_dp, s_cnt(cat) == 1)
                std_val = 0.0_dp
            end if
            write(25, '(2F12.4)', advance='no') avg_val, std_val
        end do
        write(25, *)
    end do

    ! --- メイン解析ループ2: O原子ベース (d5 / HB coordination) ---
    print *, "Analyzing Oxygen-base properties (d5, HB-coordination)..."
    do t = t_init, t_last
        current_hb_count = 0
        ! d5計算
        do i_o = ions(1,1), ions(2,1)
            block
                real(dp) :: dists(5), t_d
                integer :: m
                dists = 999.9_dp
                do j_o = ions(1,1), ions(2,1)
                    if (i_o == j_o) cycle
                    call get_dist_pbc(coords(:,i_o,t), coords(:,j_o,t), cell_data(:,t), t_d)
                    if (t_d < dists(5)) then
                        dists(5) = t_d
                        do m = 4, 1, -1
                            if (dists(m+1) < dists(m)) then
                                t_d = dists(m); dists(m) = dists(m+1); dists(m+1) = t_d
                            else; exit; end if
                        end do
                    end if
                end do
                avg_d5_dist(i_o) = avg_d5_dist(i_o) + dists(5)
            end block
        end do
        ! 酸素周りのHB本数カウント
        do i_h = ions(1,2), ions(2,2)
            d1 = 999.9_dp; d2 = 999.9_dp
            do i_o = ions(1,1), ions(2,1)
                call get_dist_pbc(coords(:,i_h,t), coords(:,i_o,t), cell_data(:,t), datoms)
                if (datoms < d1) then
                    d2 = d1; snd = fst; d1 = datoms; fst = i_o
                else if (datoms < d2) then
                    d2 = datoms; snd = i_o
                end if
            end do
            call get_angle_3pt(coords(:,i_h,t), coords(:,fst,t), coords(:,snd,t), cell_data(:,t), h_ang)
            if (d2 < 2.5_dp .and. (PI - h_ang) < (35.0_dp*PI/180.0_dp)) then
                current_hb_count(fst) = current_hb_count(fst) + 1
                current_hb_count(snd) = current_hb_count(snd) + 1
            end if
        end do
        do i_o = ions(1,1), ions(2,1)
            j = min(6, current_hb_count(i_o))
            hb_coord_distri(j) = hb_coord_distri(j) + 1.0_dp
        end do
    end do

    ! --- 最終出力 ---
    write(23, '(A5, 4X, A7)') "r (A)", "d5_rate"
    do s = 0, 300
        r_threshold = 2.0_dp + real(s, dp) * 0.01_dp
        i = count( (avg_d5_dist / real(step_count, dp)) < r_threshold )
        write(23, '(F5.2, 4X, F5.2)') r_threshold, real(i, dp)/real(n_o, dp)
    end do
    
    do j = 0, 6
        write(24, '(I2, F12.6)') j, hb_coord_distri(j) / (real(step_count, dp) * real(n_o, dp))
    end do

    close(21); close(22); close(23); close(24); close(25)
    print *, "All analysis completed successfully."

end program