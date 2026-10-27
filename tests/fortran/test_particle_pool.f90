!> 粒子 batch の生成と追加、バッチ注入 layout、shared particle pool（rank 共有 slot pool）を検証するテスト。
program test_particle_pool
  use bem_kinds, only: dp, i32, i64
  use bem_mesh, only: init_mesh
  use bem_simulator, only: run_absorption_insulator
  use bem_app_config, only: app_config, default_app_config, species_from_defaults, seed_particles_from_config, &
                            init_particle_batch_from_config, particle_inflow_reservoir, &
                            batch_injection_layout_type, build_batch_injection_layout, generate_batch_slot, &
                            slot_kind_volume, slot_kind_reservoir_face
  use bem_app_config_runtime, only: particle_source_plan_type, build_particle_source_plan
  use bem_types, only: mesh_type, particles_soa, sim_stats, injection_state, bc_open, bc_reflect
  use bem_charge_ledger, only: charge_ledger_type
  use bem_injection, only: sample_reservoir_face_particle_counter
  use bem_particles, only: init_particles, append_particles
  use bem_electrostatic_snapshot, only: electrostatic_snapshot_type
  use bem_external_boundary_contract, only: external_inflow_scalar_barrier
  use test_support, only: test_init, test_begin, test_end, test_summary, &
                          assert_true, assert_equal_i32, assert_equal_i64, assert_close_dp, assert_allclose_1d
  implicit none

  type(app_config) :: cfg
  type(particle_source_plan_type) :: source_plan
  type(batch_injection_layout_type) :: layout
  type(injection_state) :: state, state_static
  type(particles_soa) :: pcls
  real(dp) :: x(3), v(3), w, q, m, x2(3), v2(3), w2, q2, m2
  integer(i32) :: species_idx, source_element, collision_status, bounce_count
  integer(i32) :: species_idx2, source_element2, collision_status2, bounce_count2
  logical :: emitted, emitted2
  integer(i32) :: slot

  call test_init(8)

  call test_begin('layout_segments_follow_static_order')
  call setup_layout_fixture(cfg)
  allocate (state%macro_residual(2), state%boundary_macro_residual(6, 2))
  state%macro_residual = 0.0_dp
  state%boundary_macro_residual = 0.0_dp
  call build_particle_source_plan(cfg, source_plan)
  call build_batch_injection_layout(cfg, 1_i32, layout, state=state, source_plan=source_plan)
  call assert_true(layout%ready, 'layout should be ready')
  call assert_equal_i32(layout%total_slots, 6_i32, 'layout total slot count mismatch')
  call assert_equal_i32(layout%nsegments, 2_i32, 'layout segment count mismatch')
  call assert_equal_i32(layout%segment_species(1), 1_i32, 'segment 1 species mismatch')
  call assert_equal_i32(layout%segment_kind(1), slot_kind_volume, 'segment 1 kind mismatch')
  call assert_equal_i32(layout%segment_begin(1), 1_i32, 'segment 1 begin mismatch')
  call assert_equal_i32(layout%segment_end(1), 5_i32, 'segment 1 end mismatch')
  call assert_equal_i32(layout%segment_species(2), 2_i32, 'segment 2 species mismatch')
  call assert_equal_i32(layout%segment_kind(2), slot_kind_reservoir_face, 'segment 2 kind mismatch')
  call assert_equal_i32(layout%segment_face(2), 5_i32, 'segment 2 face mismatch')
  call assert_equal_i32(layout%segment_begin(2), 6_i32, 'segment 2 begin mismatch')
  call assert_equal_i32(layout%segment_end(2), 6_i32, 'segment 2 end mismatch')
  call assert_equal_i32(layout%boundary_global_count(5, 2), 1_i32, 'reservoir global count mismatch')
  ! Slots interleave the segments round-robin: volume, reservoir, volume, volume, ...
  call assert_equal_i32(layout%slot_segment(1), 1_i32, 'slot 1 should come from the volume segment')
  call assert_equal_i32(layout%slot_segment(2), 2_i32, 'slot 2 should come from the reservoir segment')
  call assert_equal_i32(layout%slot_index(2), 1_i32, 'slot 2 segment index mismatch')
  call assert_equal_i32(layout%slot_segment(6), 1_i32, 'slot 6 should come from the volume segment')
  call assert_equal_i32(layout%slot_index(6), 5_i32, 'slot 6 segment index mismatch')
  ! The static path samples the same global counts on one rank.
  allocate (state_static%macro_residual(2), state_static%boundary_macro_residual(6, 2))
  state_static%macro_residual = 0.0_dp
  state_static%boundary_macro_residual = 0.0_dp
  call seed_particles_from_config(cfg)
  call init_particle_batch_from_config(cfg, 1_i32, pcls, state=state_static, source_plan=source_plan)
  call assert_equal_i32(pcls%n, layout%total_slots, 'static batch should hold every layout slot on one rank')
  call assert_equal_i32(pcls%species_id(2), 2_i32, 'static batch round-robin should interleave species')
  call assert_close_dp( &
    state_static%boundary_macro_residual(5, 2), state%boundary_macro_residual(5, 2), 1.0e-15_dp, &
    'layout and static path must advance the reservoir residual identically' &
    )
  call test_end()

  call test_begin('slot_generation_is_deterministic')
  do slot = 1_i32, layout%total_slots
    call generate_batch_slot( &
      cfg, layout, slot, x=x, v=v, w=w, q=q, m=m, species_idx=species_idx, source_element=source_element, &
      emitted=emitted, collision_status=collision_status, bounce_count=bounce_count &
      )
    call generate_batch_slot( &
      cfg, layout, slot, x=x2, v=v2, w=w2, q=q2, m=m2, species_idx=species_idx2, source_element=source_element2, &
      emitted=emitted2, collision_status=collision_status2, bounce_count=bounce_count2 &
      )
    call assert_true(emitted .and. emitted2, 'volume and reservoir slots always emit')
    call assert_true(all(x == x2) .and. all(v == v2) .and. w == w2 .and. q == q2 .and. m == m2, &
                     'the same slot must regenerate the same particle')
    call assert_equal_i32(species_idx, species_idx2, 'the same slot must regenerate the same species')
    call assert_equal_i32(source_element, -1_i32, 'volume and reservoir slots have no source element')
    if (layout%slot_segment(slot) == 1_i32) then
      call assert_equal_i32(species_idx, 1_i32, 'volume slot species mismatch')
      call assert_true(all(x >= cfg%particle_species(1)%pos_low) .and. all(x <= cfg%particle_species(1)%pos_high), &
                       'volume slot position must lie in the seed box')
      call assert_close_dp(w, 2.0_dp, 0.0_dp, 'volume slot weight mismatch')
    else
      call assert_equal_i32(species_idx, 2_i32, 'reservoir slot species mismatch')
      call assert_close_dp(x(3), nearest(0.0_dp, 1.0_dp), 0.0_dp, 'reservoir slot must start on the z_low face')
      call assert_true(v(3) > 0.0_dp, 'reservoir slot must move inward')
      call assert_true(x(1) >= 0.0_dp .and. x(1) <= 1.0_dp .and. x(2) >= 0.0_dp .and. x(2) <= 1.0_dp, &
                       'reservoir slot must lie on the face')
    end if
  end do
  ! Different slots of one segment draw different numbers.
  call generate_batch_slot( &
    cfg, layout, 1_i32, x=x, v=v, w=w, q=q, m=m, species_idx=species_idx, source_element=source_element, &
    emitted=emitted, collision_status=collision_status, bounce_count=bounce_count &
    )
  call generate_batch_slot( &
    cfg, layout, 3_i32, x=x2, v=v2, w=w2, q=q2, m=m2, species_idx=species_idx2, source_element=source_element2, &
    emitted=emitted2, collision_status=collision_status2, bounce_count=bounce_count2 &
    )
  call assert_true(any(x /= x2), 'different slots must sample different positions')
  call test_end()

  call test_begin('counter_reservoir_sampler_barrier_shift')
  call run_reservoir_counter_sampler_test()
  call test_end()

  call test_begin('shared_pool_matches_static_simulation')
  call run_shared_pool_simulation_test()
  call test_end()

  call test_begin('shared_pool_photo_rays_independent_of_chunk')
  call run_shared_pool_photo_test()
  call test_end()

  call test_begin('particle_append_preserves_existing_state')
  call test_particle_append_preserves_existing_state()
  call test_end()

  call test_begin('uneven_species_batch_attributes_and_empty_batch')
  call test_uneven_species_batch()
  call test_end()

  call test_begin('photo_batch')
  call test_photo_batch()
  call test_end()

  call test_summary()

contains

  subroutine test_particle_append_preserves_existing_state()
    type(particles_soa) :: particles
    real(dp) :: initial_x(3, 1), initial_v(3, 1), append_x(3, 2), append_v(3, 2)

    initial_x(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp]
    initial_v(:, 1) = [4.0_dp, 5.0_dp, 6.0_dp]
    append_x(:, 1) = [7.0_dp, 8.0_dp, 9.0_dp]
    append_x(:, 2) = [10.0_dp, 11.0_dp, 12.0_dp]
    append_v(:, 1) = [13.0_dp, 14.0_dp, 15.0_dp]
    append_v(:, 2) = [16.0_dp, 17.0_dp, 18.0_dp]
    call init_particles( &
      particles, initial_x, initial_v, [1.0_dp], [2.0_dp], [3.0_dp], &
      species_id=[1_i32], source_element=[7_i32] &
      )
    particles%alive(1) = .false.
    call append_particles( &
      particles, append_x, append_v, [4.0_dp, 5.0_dp], [6.0_dp, 7.0_dp], [8.0_dp, 9.0_dp], &
      [2_i32, 3_i32], source_element=[8_i32, 9_i32] &
      )

    call assert_equal_i32(particles%n, 3_i32, 'particle append count mismatch')
    call assert_true(.not. particles%alive(1), 'particle append must preserve prior alive state')
    call assert_true(all(particles%alive(2:3)), 'particle append must activate appended particles')
    call assert_close_dp(particles%x(3, 3), 12.0_dp, 0.0_dp, 'particle append position mismatch')
    call assert_close_dp(particles%w(3), 9.0_dp, 0.0_dp, 'particle append weight mismatch')
    call assert_equal_i32(particles%species_id(3), 3_i32, 'particle append species mismatch')
    call assert_true( &
      all(particles%source_element == [7_i32, 8_i32, 9_i32]), &
      'particle append must preserve and append source provenance' &
      )
  end subroutine test_particle_append_preserves_existing_state

  subroutine test_uneven_species_batch()
    type(app_config) :: app
    type(particle_source_plan_type) :: plan
    type(particles_soa) :: particles
    integer(i32), parameter :: counts(4) = [4_i32, 20_i32, 2_i32, 1_i32]
    integer :: s, j

    call default_app_config(app)
    app%sim%batch_count = 1_i32
    app%n_particle_species = 4_i32
    do s = 1, 4
      app%particle_species(s) = species_from_defaults()
      app%particle_species(s)%npcls_per_step = counts(s)
      app%particle_species(s)%enabled = s /= 2
      app%particle_species(s)%pos_low = real(s, dp)
      app%particle_species(s)%pos_high = real(s, dp)
      app%particle_species(s)%drift_velocity = -real(s, dp)
      app%particle_species(s)%temperature_k = 0.0_dp
      app%particle_species(s)%q_particle = -real(s, dp)
      app%particle_species(s)%m_particle = 2.0_dp*s
      app%particle_species(s)%w_particle = 3.0_dp*s
    end do
    call seed_particles_from_config(app)
    call build_particle_source_plan(app, plan)
    call init_particle_batch_from_config(app, 1_i32, particles, source_plan=plan)
    ! The interleaving order is layout_segments_follow_static_order's contract; here every particle must carry
    ! the attributes of its own species, and disabled species emit nothing.
    call assert_equal_i32(particles%n, 7_i32, 'uneven species count mismatch')
    do s = 1, 4
      call assert_equal_i32(int(count(particles%species_id == s), i32), merge(0_i32, counts(s), s == 2), &
                            'per-species batch count mismatch')
    end do
    do j = 1, particles%n
      s = particles%species_id(j)
      call assert_true(all(particles%x(:, j) == real(s, dp)), 'position was assigned to the wrong species')
      call assert_true(all(particles%v(:, j) == -real(s, dp)), 'velocity was assigned to the wrong species')
      call assert_close_dp(particles%q(j), -real(s, dp), 0.0_dp, 'species charge mismatch')
      call assert_close_dp(particles%m(j), 2.0_dp*s, 0.0_dp, 'species mass mismatch')
      call assert_close_dp(particles%w(j), 3.0_dp*s, 0.0_dp, 'species weight mismatch')
    end do
    call assert_true(all(particles%alive), 'new batch particles must be alive')
    call assert_true(all(particles%source_element == -1_i32), 'volume seed acquired surface provenance')
    app%particle_species(:)%npcls_per_step = 0_i32
    call build_particle_source_plan(app, plan)
    call init_particle_batch_from_config(app, 1_i32, particles, source_plan=plan)
    call assert_equal_i32(particles%n, 0_i32, 'empty batch retained old particles')
    call assert_true(size(particles%x, 1) == 3 .and. size(particles%x, 2) == 0, 'empty batch position shape mismatch')
    call assert_true(size(particles%species_id) == 0 .and. size(particles%alive) == 0, 'empty batch state mismatch')
  end subroutine test_uneven_species_batch

  !> volume、密度 0 の reservoir、photo_raycast が混在する batch の粒子数・電荷・光電子の反作用電荷の配分。
  subroutine test_photo_batch()
    type(mesh_type) :: mesh
    type(app_config) :: cfg
    type(particle_source_plan_type) :: plan
    type(particles_soa) :: pcls
    type(injection_state) :: state
    type(electrostatic_snapshot_type) :: snapshot
    real(dp), allocatable :: photo_emission_dq(:), photo_source_charge(:)
    real(dp) :: photo_v0(3, 1), photo_v1(3, 1), photo_v2(3, 1)
    real(dp) :: expected_photo_counter, observed_photo_counter
    integer(i32) :: i

    photo_v0(:, 1) = [0.0d0, 0.0d0, 0.1d0]
    photo_v1(:, 1) = [1.0d0, 0.0d0, 0.1d0]
    photo_v2(:, 1) = [0.0d0, 1.0d0, 0.1d0]
    call init_mesh(mesh, photo_v0, photo_v1, photo_v2)
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals

    call default_app_config(cfg)
    cfg%sim%rng_seed = 999_i32
    cfg%sim%batch_count = 1_i32
    cfg%sim%batch_duration = 1.0d0
    cfg%sim%use_box = .true.
    cfg%sim%box_min = [0.0d0, 0.0d0, 0.0d0]
    cfg%sim%box_max = [1.0d0, 1.0d0, 1.0d0]
    cfg%sim%reservoir_potential_model = 'infinity_barrier'
    cfg%boundary_contract%inflow_map = external_inflow_scalar_barrier
    cfg%sim%phi_infty = 0.0d0
    cfg%sim%injection_face_phi_grid_n = 2_i32
    cfg%n_particle_species = 3_i32

    cfg%particle_species(1) = species_from_defaults()
    cfg%particle_species(1)%source_mode = 'volume_seed'
    cfg%particle_species(1)%npcls_per_step = 1_i32
    cfg%particle_species(1)%q_particle = 1.0d0
    cfg%particle_species(1)%m_particle = 1.0d0
    cfg%particle_species(1)%w_particle = 1.0d0
    cfg%particle_species(1)%pos_low = [0.0d0, 0.0d0, 0.5d0]
    cfg%particle_species(1)%pos_high = [0.0d0, 0.0d0, 0.5d0]
    cfg%particle_species(1)%drift_velocity = [0.0d0, 0.0d0, 0.0d0]
    cfg%particle_species(1)%temperature_k = 0.0d0

    cfg%particle_species(2) = species_from_defaults()
    cfg%particle_species(2)%source_mode = 'volume_seed'
    cfg%particle_species(2)%npcls_per_step = 0_i32
    cfg%particle_species(2)%number_density_cm3 = 0.0d0
    cfg%particle_species(2)%temperature_k = 0.0d0
    cfg%particle_species(2)%q_particle = 1.0d0
    cfg%particle_species(2)%m_particle = 1.0d0
    cfg%particle_species(2)%w_particle = 10.0d0
    cfg%particle_species(2)%boundary_inflow_low(3) = particle_inflow_reservoir
    cfg%particle_species(2)%drift_velocity = [0.0d0, 0.0d0, 1.0d0]

    cfg%particle_species(3) = species_from_defaults()
    cfg%particle_species(3)%source_mode = 'photo_raycast'
    cfg%particle_species(3)%emit_current_density_a_m2 = 1.0d0
    cfg%particle_species(3)%rays_per_batch = 1_i32
    cfg%particle_species(3)%normal_drift_speed = 0.0d0
    cfg%particle_species(3)%ray_direction = [0.0d0, 0.0d0, -1.0d0]
    cfg%particle_species(3)%has_ray_direction = .true.
    cfg%particle_species(3)%q_particle = -1.0d0
    cfg%particle_species(3)%m_particle = 1.0d0
    cfg%particle_species(3)%temperature_k = 0.0d0
    cfg%particle_species(3)%deposit_opposite_charge_on_emit = .true.
    cfg%particle_species(3)%inject_face = 'z_high'
    cfg%particle_species(3)%pos_low = [0.0d0, 0.0d0, 1.0d0]
    cfg%particle_species(3)%pos_high = [0.1d0, 0.1d0, 1.0d0]

    allocate (state%macro_residual(3), state%boundary_macro_residual(6, 3))
    state%boundary_macro_residual = 0.0_dp
    state%macro_residual = 0.0d0
    allocate (photo_emission_dq(mesh%nelem), photo_source_charge(mesh%nelem))
    call snapshot%init(mesh, cfg%sim, cfg%field, cfg%periodic2, cfg%panel)
    call snapshot%refresh(mesh)
    call seed_particles_from_config(cfg)
    call build_particle_source_plan(cfg, plan)
    call init_particle_batch_from_config( &
      cfg, 1_i32, pcls, state, mesh=mesh, snapshot=snapshot, photo_emission_dq=photo_emission_dq, source_plan=plan &
      )
    call assert_equal_i32(pcls%n, 2_i32, 'mixed batch particle count mismatch')
    call assert_equal_i32(int(count(pcls%species_id == 1_i32), i32), 1_i32, 'volume species count mismatch')
    call assert_equal_i32(int(count(pcls%species_id == 2_i32), i32), 0_i32, 'zero-density reservoir count mismatch')
    call assert_equal_i32(int(count(pcls%species_id == 3_i32), i32), 1_i32, 'photo species count mismatch')
    call assert_close_dp( &
      sum(pack(pcls%q, pcls%species_id == 1_i32)), 1.0d0, 1.0d-12, &
      'mixed batch species-1 charge mismatch' &
      )
    expected_photo_counter = cfg%particle_species(3)%emit_current_density_a_m2* &
                             (cfg%particle_species(3)%pos_high(1) - cfg%particle_species(3)%pos_low(1))* &
                             (cfg%particle_species(3)%pos_high(2) - cfg%particle_species(3)%pos_low(2))* &
                             cfg%sim%batch_duration
    observed_photo_counter = 0.0_dp
    photo_source_charge = 0.0_dp
    do i = 1_i32, pcls%n
      if (pcls%species_id(i) == 3_i32) then
        call assert_true(pcls%source_element(i) >= 1_i32 .and. pcls%source_element(i) <= mesh%nelem, &
                         'photo particle source provenance is outside the emitting mesh')
        observed_photo_counter = observed_photo_counter - pcls%q(i)*pcls%w(i)
        photo_source_charge(pcls%source_element(i)) = photo_source_charge(pcls%source_element(i)) - &
                                                      pcls%q(i)*pcls%w(i)
        call assert_true(pcls%w(i) > 0.0d0, 'mixed batch photo weight should be positive')
      else
        call assert_equal_i32(pcls%source_element(i), -1_i32, &
                              'non-photo batch particle unexpectedly carries surface provenance')
      end if
    end do
    call assert_close_dp(observed_photo_counter, expected_photo_counter, 1.0d-14, 'photo current budget mismatch')
    call assert_close_dp(sum(photo_emission_dq), expected_photo_counter, 1.0d-12, 'photo counter charge mismatch')
    call assert_allclose_1d(photo_source_charge, photo_emission_dq, 1.0d-12, &
                            'source provenance does not reproduce the emitted countercharge distribution')
    call assert_close_dp(state%boundary_macro_residual(5, 2), 0.0_dp, 0.0_dp, 'zero-density reservoir residual mismatch')
  end subroutine test_photo_batch

  !> species 1: 5 個の volume_seed、species 2: z_low からの reservoir 流入（1 バッチ 1 個）。
  subroutine setup_layout_fixture(cfg)
    type(app_config), intent(out) :: cfg

    call default_app_config(cfg)
    cfg%sim%rng_seed = 4242_i32
    cfg%sim%batch_count = 2_i32
    cfg%sim%batch_duration = 1.0_dp
    cfg%sim%has_batch_duration = .true.
    cfg%sim%use_box = .true.
    cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    cfg%n_particle_species = 2_i32
    cfg%particle_species(1) = species_from_defaults()
    cfg%particle_species(1)%source_mode = 'volume_seed'
    cfg%particle_species(1)%npcls_per_step = 5_i32
    cfg%particle_species(1)%q_particle = 1.0_dp
    cfg%particle_species(1)%m_particle = 1.0_dp
    cfg%particle_species(1)%w_particle = 2.0_dp
    cfg%particle_species(1)%pos_low = [0.2_dp, 0.3_dp, 0.4_dp]
    cfg%particle_species(1)%pos_high = [0.6_dp, 0.7_dp, 0.8_dp]
    cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -1.0_dp]
    cfg%particle_species(1)%temperature_k = 1.0_dp
    cfg%particle_species(2) = species_from_defaults()
    cfg%particle_species(2)%source_mode = 'volume_seed'
    cfg%particle_species(2)%npcls_per_step = 0_i32
    cfg%particle_species(2)%number_density_m3 = 1.0_dp
    cfg%particle_species(2)%temperature_k = 0.0_dp
    cfg%particle_species(2)%q_particle = 0.0_dp
    cfg%particle_species(2)%m_particle = 1.0_dp
    cfg%particle_species(2)%w_particle = 1.0_dp
    cfg%particle_species(2)%has_w_particle = .true.
    cfg%particle_species(2)%boundary_inflow_low(3) = particle_inflow_reservoir
    cfg%particle_species(2)%drift_velocity = [0.0_dp, 0.0_dp, 1.0_dp]
  end subroutine setup_layout_fixture

  !> 固定位置・零温度の 4 粒子が 1 要素に吸収される fixture（粒子は乱数に依らない）。
  subroutine setup_absorption_fixture(mesh, cfg)
    type(mesh_type), intent(out) :: mesh
    type(app_config), intent(out) :: cfg
    real(dp) :: v0(3, 1), v1(3, 1), v2(3, 1)

    v0(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    v1(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    v2(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
    call init_mesh(mesh, v0, v1, v2)
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals

    call default_app_config(cfg)
    cfg%sim%rng_seed = 2468_i32
    cfg%sim%batch_count = 2_i32
    cfg%sim%dt = 1.0_dp
    cfg%sim%max_step = 1_i32
    cfg%sim%q_floor = 1.0e-30_dp
    cfg%sim%use_box = .true.
    cfg%sim%box_min = [-1.0_dp, -1.0_dp, -2.0_dp]
    cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    cfg%n_particle_species = 1_i32
    cfg%particle_species(1) = species_from_defaults()
    cfg%particle_species(1)%source_mode = 'volume_seed'
    cfg%particle_species(1)%npcls_per_step = 4_i32
    cfg%particle_species(1)%q_particle = 1.0_dp
    cfg%particle_species(1)%m_particle = 1.0_dp
    cfg%particle_species(1)%w_particle = 1.0_dp
    cfg%particle_species(1)%pos_low = [0.2_dp, 0.2_dp, 0.8_dp]
    cfg%particle_species(1)%pos_high = [0.2_dp, 0.2_dp, 0.8_dp]
    cfg%particle_species(1)%drift_velocity = [0.0_dp, 0.0_dp, -2.0_dp]
    cfg%particle_species(1)%temperature_k = 0.0_dp
  end subroutine setup_absorption_fixture

  !> 零温度では法線速度がドリフトに固定されるので、障壁エネルギーの符号付きシフトを解析値 sqrt(u^2 - barrier) と比べられる。
  subroutine run_reservoir_counter_sampler_test()
    real(dp) :: x_slot(3), v_slot(3)
    real(dp), parameter :: box_min(3) = [0.0_dp, 0.0_dp, 0.0_dp], box_max(3) = [1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: drift(3) = [0.0_dp, 0.0_dp, 1.0_dp]
    real(dp) :: barrier
    integer :: trial

    do trial = 1, 2
      ! An accelerating sheath (negative barrier) raises the speed; a retarding one lowers it.
      barrier = merge(-3.0_dp, 0.75_dp, trial == 1)
      call sample_reservoir_face_particle_counter( &
        box_min, box_max, 'z_low', box_min, [1.0_dp, 1.0_dp, 0.0_dp], drift, 1.0_dp, 0.0_dp, barrier, 0.0_dp, &
        cfg%sim%rng_seed, [7_i64, 1_i64], x_slot, v_slot &
        )
      call assert_close_dp(v_slot(3), sqrt(1.0_dp - barrier), 1.0e-12_dp, 'counter sampler barrier shift mismatch')
      call assert_close_dp(x_slot(3), nearest(0.0_dp, 1.0_dp), 0.0_dp, 'counter sampler must start one ulp inside the face')
    end do
  end subroutine run_reservoir_counter_sampler_test

  subroutine run_shared_pool_simulation_test()
    type(mesh_type) :: mesh_static, mesh_shared
    type(app_config) :: cfg_static, cfg_shared
    type(sim_stats) :: stats_static, stats_shared
    type(charge_ledger_type) :: ledger_static, ledger_shared

    call setup_absorption_fixture(mesh_static, cfg_static)
    call seed_particles_from_config(cfg_static)
    call run_absorption_insulator(mesh_static, cfg_static, stats_static, charge_ledger=ledger_static)
    call assert_equal_i64(stats_static%processed_particles, 8_i64, 'static processed_particles mismatch')
    call assert_equal_i64(stats_static%particle_pool_claims, 0_i64, 'static mode must not claim pool slots')

    call setup_absorption_fixture(mesh_shared, cfg_shared)
    cfg_shared%sim%particle_distribution = 'shared'
    cfg_shared%sim%particle_pool_chunk = 3_i32
    call seed_particles_from_config(cfg_shared)
    call run_absorption_insulator(mesh_shared, cfg_shared, stats_shared, charge_ledger=ledger_shared)
    call assert_equal_i64(stats_shared%processed_particles, 8_i64, 'shared processed_particles mismatch')
    call assert_equal_i64(stats_shared%absorbed, 8_i64, 'shared absorbed mismatch')
    call assert_equal_i64(stats_shared%escaped, 0_i64, 'shared escaped mismatch')
    call assert_equal_i32(stats_shared%batches, 2_i32, 'shared batches mismatch')
    call assert_close_dp(mesh_shared%q_elem(1), mesh_static%q_elem(1), 1.0e-12_dp, 'shared deposited charge mismatch')
    call assert_equal_i64(ledger_shared%injected_count(1), 8_i64, 'shared ledger injected count mismatch')
    call assert_equal_i64(ledger_shared%absorbed_count(1), 8_i64, 'shared ledger absorbed count mismatch')
    call assert_close_dp(ledger_shared%residual(), 0.0_dp, 1.0e-12_dp, 'shared ledger residual mismatch')
    ! 4 slots in chunks of 3 take at least 2 claims per batch plus one empty claim per thread.
    call assert_true(stats_shared%particle_pool_claims >= 4_i64, 'shared mode must record its pool claims')
    call assert_true(stats_shared%particle_pool_claim_seconds >= 0.0_dp, 'claim seconds must be nonnegative')
  end subroutine run_shared_pool_simulation_test

  !> photo ray の pool では、外れた ray が slot を空けるので batch は claim した分より小さく詰め直される。
  !! chunk 幅を変えても同じ ray 集合が同じ要素に当たるので、収支は丸め誤差の範囲で一致する。
  subroutine run_shared_pool_photo_test()
    integer(i32), parameter :: rays = 1024_i32
    type(mesh_type) :: mesh_small, mesh_large
    type(app_config) :: cfg_small, cfg_large
    type(sim_stats) :: stats_small, stats_large
    type(charge_ledger_type) :: ledger_small, ledger_large
    real(dp) :: signature_small(9), signature_large(9)
    integer :: value_index
    character(len=96) :: message

    call setup_photo_fixture(mesh_small, cfg_small, rays)
    cfg_small%sim%particle_distribution = 'shared'
    cfg_small%sim%particle_pool_chunk = 7_i32
    call seed_particles_from_config(cfg_small)
    call run_absorption_insulator(mesh_small, cfg_small, stats_small, charge_ledger=ledger_small)

    call setup_photo_fixture(mesh_large, cfg_large, rays)
    cfg_large%sim%particle_distribution = 'shared'
    cfg_large%sim%particle_pool_chunk = 4096_i32
    call seed_particles_from_config(cfg_large)
    call run_absorption_insulator(mesh_large, cfg_large, stats_large, charge_ledger=ledger_large)

    call assert_true(stats_small%processed_particles > 0_i64, 'photo pool must emit particles')
    call assert_true(stats_small%processed_particles < int(rays, i64), 'rays that miss must not become particles')
    call assert_equal_i64(stats_large%processed_particles, stats_small%processed_particles, &
                          'photo pool particle count must not depend on the chunk width')
    call assert_equal_i64(ledger_large%emitted_count(1), ledger_small%emitted_count(1), &
                          'photo pool emitted count must not depend on the chunk width')
    call pack_photo_signature(mesh_small, ledger_small, signature_small)
    call pack_photo_signature(mesh_large, ledger_large, signature_large)
    do value_index = 1, size(signature_small)
      write (message, '(a,i0)') 'photo pool signature mismatch at index ', value_index
      call assert_close_dp(signature_large(value_index), signature_small(value_index), &
                           1.0e-11_dp*max(1.0_dp, abs(signature_small(value_index))), trim(message))
    end do
  end subroutine run_shared_pool_photo_test

  subroutine setup_photo_fixture(mesh, cfg, rays_per_batch)
    type(mesh_type), intent(out) :: mesh
    type(app_config), intent(out) :: cfg
    integer(i32), intent(in) :: rays_per_batch
    real(dp) :: tri_v0(3, 4), tri_v1(3, 4), tri_v2(3, 4)

    ! Two squares (four triangles) cover half of the z=0.2 plane, so half of the rays miss and escape.
    tri_v0(:, 1) = [0.0_dp, 0.0_dp, 0.2_dp]
    tri_v1(:, 1) = [0.5_dp, 0.0_dp, 0.2_dp]
    tri_v2(:, 1) = [0.5_dp, 0.5_dp, 0.2_dp]
    tri_v0(:, 2) = [0.0_dp, 0.0_dp, 0.2_dp]
    tri_v1(:, 2) = [0.5_dp, 0.5_dp, 0.2_dp]
    tri_v2(:, 2) = [0.0_dp, 0.5_dp, 0.2_dp]
    tri_v0(:, 3) = [0.5_dp, 0.5_dp, 0.2_dp]
    tri_v1(:, 3) = [1.0_dp, 0.5_dp, 0.2_dp]
    tri_v2(:, 3) = [1.0_dp, 1.0_dp, 0.2_dp]
    tri_v0(:, 4) = [0.5_dp, 0.5_dp, 0.2_dp]
    tri_v1(:, 4) = [1.0_dp, 1.0_dp, 0.2_dp]
    tri_v2(:, 4) = [0.5_dp, 1.0_dp, 0.2_dp]
    call init_mesh(mesh, tri_v0, tri_v1, tri_v2)
    mesh%elem_vacuum_sign = 1_i32
    mesh%vacuum_normals = mesh%normals

    call default_app_config(cfg)
    cfg%sim%rng_seed = 998_i32
    cfg%sim%batch_count = 1_i32
    cfg%sim%batch_duration = 1.0_dp
    cfg%sim%has_batch_duration = .true.
    cfg%sim%dt = 0.6_dp
    cfg%sim%max_step = 1_i32
    cfg%sim%q_floor = 1.0e-30_dp
    cfg%sim%field_solver = 'direct'
    cfg%field%backend = 'direct'
    cfg%panel%kernel_id = 'triangle_p0_exact_direct'
    cfg%sim%use_box = .true.
    cfg%sim%box_min = [0.0_dp, 0.0_dp, 0.0_dp]
    cfg%sim%box_max = [1.0_dp, 1.0_dp, 1.0_dp]
    cfg%sim%bc_high(3) = bc_open
    cfg%n_particle_species = 1_i32
    cfg%particle_species(1) = species_from_defaults()
    cfg%particle_species(1)%source_mode = 'photo_raycast'
    cfg%particle_species(1)%rays_per_batch = rays_per_batch
    cfg%particle_species(1)%emit_current_density_a_m2 = 1.0_dp
    cfg%particle_species(1)%deposit_opposite_charge_on_emit = .true.
    cfg%particle_species(1)%boundary_high(3) = bc_reflect
    cfg%particle_species(1)%q_particle = -1.0_dp
    cfg%particle_species(1)%m_particle = 1.0_dp
    cfg%particle_species(1)%temperature_k = 0.0_dp
    cfg%particle_species(1)%normal_drift_speed = 1.0_dp
    cfg%particle_species(1)%inject_face = 'z_high'
    cfg%particle_species(1)%pos_low = [0.0_dp, 0.0_dp, 1.0_dp]
    cfg%particle_species(1)%pos_high = [1.0_dp, 1.0_dp, 1.0_dp]
    cfg%particle_species(1)%ray_direction = [0.0_dp, 0.0_dp, -1.0_dp]
    cfg%particle_species(1)%has_ray_direction = .true.
  end subroutine setup_photo_fixture

  subroutine pack_photo_signature(mesh, ledger, signature)
    type(mesh_type), intent(in) :: mesh
    type(charge_ledger_type), intent(in) :: ledger
    real(dp), intent(out) :: signature(9)

    signature = [ &
                mesh%q_elem, ledger%surface_charge_before, ledger%surface_charge_after, &
                ledger%emitted_from_surface(1), ledger%absorbed_on_surface(1), ledger%residual() &
                ]
  end subroutine pack_photo_signature

end program test_particle_pool
