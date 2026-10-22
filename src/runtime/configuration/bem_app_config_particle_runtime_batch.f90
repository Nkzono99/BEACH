!> バッチ注入 layout の global slot を rank に連続区間で割り振り、この rank の粒子を生成する。
!> 各 slot は (batch, species, kind, face, index) の counter-based 乱数で生成されるので、
!> rank 数によらずバッチ全体の粒子集合は同じになる。
submodule(bem_app_config_particle_runtime) bem_app_config_particle_runtime_batch
  use bem_particles, only: allocate_particles
  implicit none
contains

  !> 指定バッチ番号に対応する、この rank の粒子バッチを生成する。
  !! @param[in] cfg 粒子種とシミュレーション条件を含むアプリ設定。
  !! @param[in] batch_idx 生成対象のバッチ番号（1始まり）。
  !! @param[out] pcls 生成したバッチ粒子群（slot 順）。
  !! @param[inout] state 境界 reservoir 流入の残差状態（必要時のみ）。
  !! @param[in] mesh 現在バッチ開始時点の電荷分布メッシュ（photo_raycast と電位補正で必要）。
  !! @param[out] photo_emission_dq photo_raycast 放出起因の要素電荷差分 `photo_emission_dq(nelem)`（省略可）。
  !! @param[out] collision_failure_status 不完全な photo collision query の status（省略時は停止）。
  !! @param[out] collision_failure_species 不完全な照会を返した最小 slot の species index。
  !! @param[out] collision_failure_ray 不完全な照会を返した ray の species 内番号。
  !! @param[out] collision_failure_bounce 不完全な照会を返した bounce index。
  !! @param[inout] snapshot refresh 済み静電 snapshot（注入電位補正の使用時に必要）。
  !! @param[in] source_plan run中に再利用する粒子source導出値（省略時は呼出し内で構築）。
  !! @param[out] photo_emission_dq_by_species species別のphoto放出反作用電荷 `photo_emission_dq_by_species(nelem, nspecies)`（省略可）。
  !! @param[out] slot_ids 各粒子の global slot 番号（省略可）。
  module procedure init_particle_batch_from_config

  integer(i32) :: local_rank, n_ranks, local_count, slot_begin, local_slot, n_emit, out_idx, failed_slot
  integer(i32) :: s, elem
  type(batch_injection_layout_type) :: layout
  real(dp), allocatable :: gx(:, :), gv(:, :), gw(:), gq(:), gm(:)
  integer(i32), allocatable :: gspecies(:), gelement(:), gstatus(:), gbounce(:)
  logical, allocatable :: gemitted(:)

  if (present(collision_failure_status)) collision_failure_status = collision_query_ok
  if (present(collision_failure_species)) collision_failure_species = huge(0_i32)
  if (present(collision_failure_ray)) collision_failure_ray = huge(0_i32)
  if (present(collision_failure_bounce)) collision_failure_bounce = 0_i32
  call resolve_parallel_rank_size(local_rank, n_ranks, mpi_rank, mpi_size, mpi, 'init_particle_batch_from_config')
  if (present(photo_emission_dq)) then
    if (.not. present(mesh)) error stop 'photo_emission_dq requires mesh in init_particle_batch_from_config.'
    if (size(photo_emission_dq) /= mesh%nelem) error stop 'photo_emission_dq size mismatch.'
    photo_emission_dq = 0.0d0
  end if
  if (present(photo_emission_dq_by_species)) then
    if (.not. present(mesh)) then
      error stop 'photo_emission_dq_by_species requires mesh in init_particle_batch_from_config.'
    end if
    if (size(photo_emission_dq_by_species, 1) /= mesh%nelem .or. &
        size(photo_emission_dq_by_species, 2) /= cfg%n_particle_species) then
      error stop 'photo_emission_dq_by_species size mismatch.'
    end if
    photo_emission_dq_by_species = 0.0_dp
  end if

  call build_batch_injection_layout( &
    cfg, batch_idx, layout, state=state, mesh=mesh, mpi_rank=mpi_rank, mpi_size=mpi_size, mpi=mpi, &
    snapshot=snapshot, source_plan=source_plan &
    )
  ! Contiguous slot range of this rank, in the same split as mpi_split_count.
  local_count = mpi_split_count(layout%total_slots, local_rank, n_ranks)
  slot_begin = local_rank*(layout%total_slots/n_ranks) + min(local_rank, modulo(layout%total_slots, n_ranks)) + 1_i32

  allocate (gx(3, local_count), gv(3, local_count), gw(local_count), gq(local_count), gm(local_count))
  allocate (gspecies(local_count), gelement(local_count), gstatus(local_count), gbounce(local_count))
  allocate (gemitted(local_count))
  !$omp parallel do default(none) schedule(dynamic, 16) &
  !$omp shared(cfg, layout, mesh, local_count, slot_begin, gx, gv, gw, gq, gm, gspecies, gelement, gemitted) &
  !$omp shared(gstatus, gbounce) private(local_slot)
  do local_slot = 1_i32, local_count
    call generate_batch_slot( &
      cfg, layout, slot_begin + local_slot - 1_i32, mesh, gx(:, local_slot), gv(:, local_slot), gw(local_slot), &
      gq(local_slot), gm(local_slot), gspecies(local_slot), gelement(local_slot), gemitted(local_slot), &
      gstatus(local_slot), gbounce(local_slot) &
      )
  end do
  !$omp end parallel do

  failed_slot = 0_i32
  do local_slot = 1_i32, local_count
    if (gstatus(local_slot) == collision_query_ok) cycle
    failed_slot = local_slot
    exit
  end do
  if (failed_slot > 0_i32) then
    call finalize_particle_batch_collision_query( &
      gstatus(failed_slot), batch_idx, gspecies(failed_slot), &
      layout%slot_index(slot_begin + failed_slot - 1_i32), gbounce(failed_slot), &
      collision_failure_status, collision_failure_species, collision_failure_ray, collision_failure_bounce &
      )
    call allocate_particles(pcls, 0_i32)
    if (present(slot_ids)) allocate (slot_ids(0))
    return
  end if

  n_emit = count(gemitted)
  call allocate_particles(pcls, n_emit)
  if (present(slot_ids)) allocate (slot_ids(n_emit))
  out_idx = 0_i32
  do local_slot = 1_i32, local_count
    if (.not. gemitted(local_slot)) cycle
    out_idx = out_idx + 1_i32
    pcls%x(:, out_idx) = gx(:, local_slot)
    pcls%v(:, out_idx) = gv(:, local_slot)
    pcls%w(out_idx) = gw(local_slot)
    pcls%q(out_idx) = gq(local_slot)
    pcls%m(out_idx) = gm(local_slot)
    pcls%species_id(out_idx) = gspecies(local_slot)
    pcls%source_element(out_idx) = gelement(local_slot)
    if (present(slot_ids)) slot_ids(out_idx) = slot_begin + local_slot - 1_i32
    s = gspecies(local_slot)
    elem = gelement(local_slot)
    if (elem < 1_i32) cycle
    if (.not. cfg%particle_species(s)%deposit_opposite_charge_on_emit) cycle
    if (.not. present(photo_emission_dq) .and. .not. present(photo_emission_dq_by_species)) cycle
    if (elem > mesh%nelem) error stop 'photo_raycast emitted invalid elem_idx.'
    if (present(photo_emission_dq)) then
      photo_emission_dq(elem) = photo_emission_dq(elem) - gq(local_slot)*gw(local_slot)
    end if
    if (present(photo_emission_dq_by_species)) then
      photo_emission_dq_by_species(elem, s) = photo_emission_dq_by_species(elem, s) - gq(local_slot)*gw(local_slot)
    end if
  end do
  end procedure init_particle_batch_from_config

  !> batch injection の不完全な photo collision query を返し、status 未要求なら serial に停止する。
  subroutine finalize_particle_batch_collision_query( &
    query_status, batch_idx, species_idx, query_ray, query_bounce, &
    collision_failure_status, collision_failure_species, collision_failure_ray, collision_failure_bounce &
    )
    integer(i32), intent(in) :: query_status, batch_idx, species_idx, query_ray, query_bounce
    integer(i32), intent(out), optional :: collision_failure_status, collision_failure_species
    integer(i32), intent(out), optional :: collision_failure_ray, collision_failure_bounce
    character(len=256) :: failure_message

    if (present(collision_failure_status)) collision_failure_status = query_status
    if (present(collision_failure_species)) collision_failure_species = species_idx
    if (present(collision_failure_ray)) collision_failure_ray = query_ray
    if (present(collision_failure_bounce)) collision_failure_bounce = query_bounce
    if (present(collision_failure_status)) return

    write (failure_message, '(a,i0,a,i0,a,i0,a,i0,a,a,a,i0)') &
      'photo_raycast collision query incomplete during batch preparation: batch=', batch_idx, &
      ' species=', species_idx, ' ray=', query_ray, ' bounce=', query_bounce, &
      ' status=', trim(particle_status_name(query_status)), ' code=', query_status
    write (error_unit, '(a)') trim(failure_message)
    flush (error_unit)
    error stop 1
  end subroutine finalize_particle_batch_collision_query

end submodule bem_app_config_particle_runtime_batch
