!> 有限画像殻の面平均の滑らかな残り s(h) を、直接積分・対称性・表の補間で確かめる。
program test_periodic_shell_mean
  use bem_kinds, only: dp, i32
  use bem_periodic_shell_mean, only: periodic_shell_mean_plan_type, periodic_shell_mean_state_type, &
                                     periodic_shell_mean_excess, periodic_shell_mean_point_excess, &
                                     build_periodic_shell_mean_plan, refresh_periodic_shell_mean_state, &
                                     eval_periodic_shell_mean
  use test_support, only: test_init, test_begin, test_end, test_summary, assert_close_dp, assert_true
  implicit none

  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: cell = 1.0_dp
  real(dp), parameter :: direct_heights(3) = [0.25_dp, 1.0_dp, 2.7_dp]
  real(dp), parameter :: smooth_heights(4) = [0.0_dp, 1.0e-7_dp, 0.05_dp, 0.8_dp]
  real(dp) :: s, dsdh, s_plus, s_minus, d_plus, d_minus, reference, h, step
  real(dp) :: source(3), vertices(3, 3, 3), charge(3), value, slope, direct, direct_slope, z
  type(periodic_shell_mean_plan_type) :: plan
  type(periodic_shell_mean_state_type) :: state
  integer :: i, k

  call test_init(4)

  call test_begin('excess_matches_direct_cell_average_of_the_image_shell')
  source = [0.3_dp, 0.62_dp, 0.0_dp]
  do i = 1, 3
    h = direct_heights(i)
    call periodic_shell_mean_excess( &
      h, -1.0_dp - source(1), 2.0_dp - source(1), -1.0_dp - source(2), 2.0_dp - source(2), cell*cell, 0.0_dp, s, dsdh &
      )
    reference = direct_shell_mean(source, h) + 2.0_dp*pi*abs(h)/(cell*cell)
    call assert_close_dp(s, reference, 1.0e-9_dp*max(1.0_dp, abs(reference)), 'shell-mean excess mismatch')
  end do
  call test_end()

  call test_begin('excess_is_even_and_smooth_through_the_source_plane')
  do i = 1, 4
    h = smooth_heights(i)
    call periodic_shell_mean_excess(h, -1.2_dp, 1.8_dp, -1.45_dp, 1.55_dp, 1.0_dp, 0.0_dp, s_plus, d_plus)
    call periodic_shell_mean_excess(-h, -1.2_dp, 1.8_dp, -1.45_dp, 1.55_dp, 1.0_dp, 0.0_dp, s_minus, d_minus)
    call assert_close_dp(s_plus, s_minus, 1.0e-13_dp, 'excess must be even in h')
    call assert_close_dp(d_plus, -d_minus, 1.0e-13_dp, 'excess slope must be odd in h')
    step = 1.0e-5_dp
    call periodic_shell_mean_excess(h + step, -1.2_dp, 1.8_dp, -1.45_dp, 1.55_dp, 1.0_dp, 0.0_dp, s_plus, d_minus)
    call periodic_shell_mean_excess(h - step, -1.2_dp, 1.8_dp, -1.45_dp, 1.55_dp, 1.0_dp, 0.0_dp, s_minus, d_minus)
    call assert_close_dp(d_plus, (s_plus - s_minus)/(2.0_dp*step), 1.0e-8_dp, 'excess slope must be the derivative')
  end do
  call periodic_shell_mean_excess(0.0_dp, -1.2_dp, 1.8_dp, -1.45_dp, 1.55_dp, 1.0_dp, 0.0_dp, s, dsdh)
  call assert_close_dp(dsdh, 0.0_dp, 0.0_dp, 'excess slope must vanish on the source plane')
  call test_end()

  call test_begin('excess_is_scale_invariant_in_raw_units')
  ! Raw 1/r units: lengths scaled by lambda scale the excess by 1/lambda and its slope by 1/lambda^2.
  call periodic_shell_mean_excess(0.4_dp, -1.3_dp, 1.7_dp, -1.1_dp, 1.9_dp, 1.0_dp, 0.0_dp, s, dsdh)
  call periodic_shell_mean_excess( &
    0.4e-4_dp, -1.3e-4_dp, 1.7e-4_dp, -1.1e-4_dp, 1.9e-4_dp, 1.0e-8_dp, 0.0_dp, s_plus, d_plus &
    )
  call assert_close_dp(s_plus*1.0e-4_dp, s, 1.0e-12_dp*abs(s), 'excess must scale as inverse length')
  call assert_close_dp(d_plus*1.0e-8_dp, dsdh, 1.0e-12_dp*abs(dsdh), 'excess slope must scale as inverse length^2')
  call test_end()

  call test_begin('height_table_reproduces_the_direct_sum')
  ! One point source and two inclined panels; the Hermite field must be the derivative of its potential.
  vertices(:, :, 1) = spread([0.2_dp, 0.3_dp, 0.4_dp], 2, 3)
  vertices(:, 1, 2) = [0.1_dp, 0.1_dp, 0.0_dp]
  vertices(:, 2, 2) = [0.6_dp, 0.2_dp, 0.3_dp]
  vertices(:, 3, 2) = [0.3_dp, 0.7_dp, 0.1_dp]
  vertices(:, 1, 3) = [0.8_dp, 0.8_dp, 1.2_dp]
  vertices(:, 2, 3) = [0.9_dp, 0.4_dp, 1.0_dp]
  vertices(:, 3, 3) = [0.5_dp, 0.9_dp, 1.4_dp]
  charge = [1.0_dp, -0.7_dp, 0.4_dp]
  call build_periodic_shell_mean_plan( &
    [1_i32, 2_i32], 3_i32, [0.5_dp, 0.5_dp], [cell, cell], 1_i32, 0.0_dp, 0.0_dp, 3.0_dp, 1.0_dp/32.0_dp, &
    vertices, plan &
    )
  call refresh_periodic_shell_mean_state(plan, charge, state)
  do k = 1, 7
    z = 0.05_dp + 0.41_dp*real(k, dp)
    call eval_periodic_shell_mean(plan, state, z, value, slope)
    direct = 0.0_dp
    direct_slope = 0.0_dp
    do i = 1, 3
      call panel_excess(vertices(:, :, i), z, s, dsdh)
      direct = direct + charge(i)*s
      direct_slope = direct_slope + charge(i)*dsdh
    end do
    call assert_close_dp(value, direct, 1.0e-7_dp, 'tabulated excess mismatch')
    call assert_close_dp(slope, direct_slope, 1.0e-5_dp, 'tabulated excess slope mismatch')
    step = 1.0e-6_dp
    call eval_periodic_shell_mean(plan, state, z + step, s_plus, d_plus)
    call eval_periodic_shell_mean(plan, state, z - step, s_minus, d_minus)
    call assert_close_dp(slope, (s_plus - s_minus)/(2.0_dp*step), 1.0e-7_dp, 'tabulated slope must be the derivative')
  end do
  call assert_true(plan%nz == 97_i32, 'height table must resolve L/32')
  call test_end()

  call test_summary()

contains

  !> 3x3 画像殻の電位をセル [0,1]^2 で直接平均する（各 0.05 の正方形に 8x8 Gauss-Legendre）。
  real(dp) function direct_shell_mean(src, height) result(mean)
    real(dp), intent(in) :: src(3), height
    real(dp), parameter :: gx(8) = [ &
                           -0.9602898564975363_dp, -0.7966664774136267_dp, -0.5255324099163290_dp, -0.1834346424956498_dp, &
                           0.1834346424956498_dp, 0.5255324099163290_dp, 0.7966664774136267_dp, 0.9602898564975363_dp]
    real(dp), parameter :: gw(8) = [ &
                           0.1012285362903763_dp, 0.2223810344533745_dp, 0.3137066458778873_dp, 0.3626837833783620_dp, &
                           0.3626837833783620_dp, 0.3137066458778873_dp, 0.2223810344533745_dp, 0.1012285362903763_dp]
    integer :: a, b, p, q, ix, iy
    real(dp) :: x, y, width

    width = 0.05_dp
    mean = 0.0_dp
    do a = 0, 19
      do b = 0, 19
        do p = 1, 8
          do q = 1, 8
            x = (real(a, dp) + 0.5_dp*(gx(p) + 1.0_dp))*width
            y = (real(b, dp) + 0.5_dp*(gx(q) + 1.0_dp))*width
            do ix = -1, 1
              do iy = -1, 1
                mean = mean + 0.25_dp*width*width*gw(p)*gw(q)/ &
                       sqrt((x - src(1) - ix)**2 + (y - src(2) - iy)**2 + height**2)
              end do
            end do
          end do
        end do
      end do
    end do
  end function direct_shell_mean

  !> 表と同じ 7 点則で三角形（または点）の excess を平均する。
  subroutine panel_excess(tri, target_z, s_avg, slope_avg)
    real(dp), intent(in) :: tri(3, 3), target_z
    real(dp), intent(out) :: s_avg, slope_avg
    real(dp), parameter :: a(7) = [ &
                           1.0_dp/3.0_dp, 0.059715871789770_dp, 0.470142064105115_dp, 0.470142064105115_dp, &
                           0.797426985353087_dp, 0.101286507323456_dp, 0.101286507323456_dp]
    real(dp), parameter :: b(7) = [ &
                           1.0_dp/3.0_dp, 0.470142064105115_dp, 0.059715871789770_dp, 0.470142064105115_dp, &
                           0.101286507323456_dp, 0.797426985353087_dp, 0.101286507323456_dp]
    real(dp), parameter :: w(7) = [ &
                           0.225_dp, 0.132394152788506_dp, 0.132394152788506_dp, 0.132394152788506_dp, &
                           0.125939180544827_dp, 0.125939180544827_dp, 0.125939180544827_dp]
    real(dp) :: point(3), s_q, slope_q
    integer :: q

    s_avg = 0.0_dp
    slope_avg = 0.0_dp
    if (all(tri(:, 2) == tri(:, 1)) .and. all(tri(:, 3) == tri(:, 1))) then
      call periodic_shell_mean_point_excess(plan, tri(:, 1), target_z, s_avg, slope_avg)
      return
    end if
    do q = 1, 7
      point = a(q)*tri(:, 1) + b(q)*tri(:, 2) + (1.0_dp - a(q) - b(q))*tri(:, 3)
      call periodic_shell_mean_point_excess(plan, point, target_z, s_q, slope_q)
      s_avg = s_avg + w(q)*s_q
      slope_avg = slope_avg + w(q)*slope_q
    end do
  end subroutine panel_excess

end program test_periodic_shell_mean
