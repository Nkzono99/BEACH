module bem_outer_event_queue
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bem_kinds, only: dp, i32, i64
  use bem_interface_types, only: interface_outcome_returned_local, interface_outcome_escaped_to_infinity
  implicit none
  private

  integer(i32), parameter, public :: outer_event_outcome_return = interface_outcome_returned_local
  integer(i32), parameter, public :: outer_event_outcome_escape = interface_outcome_escaped_to_infinity

  integer(i64), parameter :: fingerprint_modulus = 2147483647_i64
  integer(i64), parameter :: fingerprint_multiplier_a = 65599_i64
  integer(i64), parameter :: fingerprint_multiplier_b = 131071_i64

  type :: fingerprint_state_type
    integer(i64) :: a = 146959810_i64
    integer(i64) :: b = 109951162_i64
  end type fingerprint_state_type

  !> An event scheduled while a macro-particle is represented by the outer domain.
  !! The record has no pointer or allocatable component so snapshots can be serialized
  !! as a flat table.  q is the physical-particle charge and w is the macro weight.
  type, public :: outer_event_record_type
    integer(i64) :: event_id = 0_i64
    real(dp) :: queued_time = 0.0_dp
    real(dp) :: due_time = 0.0_dp
    integer(i32) :: outcome = 0_i32
    integer(i32) :: species_id = 0_i32
    integer(i32) :: origin_rank = 0_i32
    integer(i32) :: origin_batch = 0_i32
    integer(i64) :: origin_particle_id = 0_i64
    integer(i32) :: interface_face_index = 0_i32
    real(dp) :: q = 0.0_dp
    real(dp) :: m = 0.0_dp
    real(dp) :: w = 0.0_dp
    real(dp) :: x(3) = 0.0_dp
    real(dp) :: v(3) = 0.0_dp
  end type outer_event_record_type

  !> Deterministic delayed-event queue ordered by (due_time, event_id).
  type, public :: outer_event_queue_type
    private
    type(outer_event_record_type), allocatable :: records(:)
    integer(i32) :: count = 0_i32
    integer(i64) :: next_event_id = 1_i64
  contains
    procedure :: init => init_outer_event_queue
    procedure :: clear => clear_outer_event_queue
    procedure :: enqueue => enqueue_outer_event
    procedure :: pop_due => pop_due_outer_events
    procedure :: size => outer_event_queue_size
    procedure :: next_id => outer_event_queue_next_id
    procedure :: snapshot => snapshot_outer_event_queue
    procedure :: restore => restore_outer_event_queue
    procedure :: signed_charge => outer_event_signed_charge
    procedure :: species_particle_number => outer_event_species_particle_number
    procedure :: photoelectron_number => outer_event_photoelectron_number
    procedure :: fingerprint_components => outer_event_queue_fingerprint_components
    procedure :: fingerprint => outer_event_queue_fingerprint
  end type outer_event_queue_type

  public :: outer_event_queue_global_fingerprint
  public :: outer_event_queue_fingerprint_is_valid

contains

  subroutine init_outer_event_queue(self, capacity, first_event_id)
    class(outer_event_queue_type), intent(out) :: self
    integer(i32), intent(in), optional :: capacity
    integer(i64), intent(in), optional :: first_event_id
    integer(i32) :: initial_capacity

    initial_capacity = 0_i32
    if (present(capacity)) initial_capacity = max(0_i32, capacity)
    self%next_event_id = 1_i64
    if (present(first_event_id)) self%next_event_id = first_event_id
    if (self%next_event_id < 1_i64) error stop 'Outer-event first ID must be positive.'
    allocate (self%records(initial_capacity))
    self%count = 0_i32
  end subroutine init_outer_event_queue

  subroutine clear_outer_event_queue(self, reset_sequence)
    class(outer_event_queue_type), intent(inout) :: self
    logical, intent(in), optional :: reset_sequence
    logical :: reset_ids

    reset_ids = .false.
    if (present(reset_sequence)) reset_ids = reset_sequence
    if (allocated(self%records) .and. self%count > 0_i32) then
      self%records(1:self%count) = outer_event_record_type()
    end if
    self%count = 0_i32
    if (reset_ids) self%next_event_id = 1_i64
  end subroutine clear_outer_event_queue

  subroutine enqueue_outer_event(self, event, event_id)
    class(outer_event_queue_type), intent(inout) :: self
    type(outer_event_record_type), intent(in) :: event
    integer(i64), intent(out), optional :: event_id
    type(outer_event_record_type) :: scheduled
    integer(i32) :: insert_at

    call validate_event_payload(event)
    if (.not. allocated(self%records)) call self%init()
    if (self%next_event_id < 1_i64 .or. self%next_event_id == huge(0_i64)) then
      error stop 'Outer-event ID sequence is exhausted.'
    end if

    scheduled = event
    scheduled%event_id = self%next_event_id
    self%next_event_id = self%next_event_id + 1_i64
    if (present(event_id)) event_id = scheduled%event_id

    call ensure_outer_event_capacity(self, self%count + 1_i32)
    insert_at = self%count + 1_i32
    do while (insert_at > 1_i32)
      if (.not. event_precedes(scheduled, self%records(insert_at - 1_i32))) exit
      self%records(insert_at) = self%records(insert_at - 1_i32)
      insert_at = insert_at - 1_i32
    end do
    self%records(insert_at) = scheduled
    self%count = self%count + 1_i32
  end subroutine enqueue_outer_event

  subroutine pop_due_outer_events(self, current_time, events)
    class(outer_event_queue_type), intent(inout) :: self
    real(dp), intent(in) :: current_time
    type(outer_event_record_type), allocatable, intent(out) :: events(:)
    integer(i32) :: due_count, old_count, remaining

    if (.not. ieee_is_finite(current_time)) error stop 'Outer-event pop time must be finite.'
    due_count = 0_i32
    do while (due_count < self%count)
      if (self%records(due_count + 1_i32)%due_time > current_time) exit
      due_count = due_count + 1_i32
    end do

    allocate (events(due_count))
    if (due_count == 0_i32) return
    events = self%records(1:due_count)

    old_count = self%count
    remaining = old_count - due_count
    if (remaining > 0_i32) self%records(1:remaining) = self%records(due_count + 1_i32:old_count)
    self%records(remaining + 1_i32:old_count) = outer_event_record_type()
    self%count = remaining
  end subroutine pop_due_outer_events

  pure integer(i32) function outer_event_queue_size(self) result(queue_size)
    class(outer_event_queue_type), intent(in) :: self

    queue_size = self%count
  end function outer_event_queue_size

  pure integer(i64) function outer_event_queue_next_id(self) result(next_id)
    class(outer_event_queue_type), intent(in) :: self

    next_id = self%next_event_id
  end function outer_event_queue_next_id

  subroutine snapshot_outer_event_queue(self, events, next_event_id)
    class(outer_event_queue_type), intent(in) :: self
    type(outer_event_record_type), allocatable, intent(out) :: events(:)
    integer(i64), intent(out), optional :: next_event_id

    allocate (events(self%count))
    if (self%count > 0_i32) events = self%records(1:self%count)
    if (present(next_event_id)) next_event_id = self%next_event_id
  end subroutine snapshot_outer_event_queue

  subroutine restore_outer_event_queue(self, events, next_event_id, capacity)
    class(outer_event_queue_type), intent(out) :: self
    type(outer_event_record_type), intent(in) :: events(:)
    integer(i64), intent(in), optional :: next_event_id
    integer(i32), intent(in), optional :: capacity
    integer(i32) :: restored_capacity, idx
    integer(i64) :: maximum_id

    if (size(events, kind=i64) > int(huge(0_i32), i64)) then
      error stop 'Outer-event snapshot is too large for the queue index kind.'
    end if
    restored_capacity = int(size(events), i32)
    if (present(capacity)) restored_capacity = max(restored_capacity, max(0_i32, capacity))
    allocate (self%records(restored_capacity))
    self%count = int(size(events), i32)
    maximum_id = 0_i64
    do idx = 1, self%count
      call validate_event_payload(events(idx))
      if (events(idx)%event_id < 1_i64) error stop 'Restored outer-event ID must be positive.'
      self%records(idx) = events(idx)
      maximum_id = max(maximum_id, events(idx)%event_id)
    end do

    call assert_unique_event_ids(self%records, self%count)
    call sort_events(self%records, self%count)

    if (present(next_event_id)) then
      self%next_event_id = next_event_id
    else if (maximum_id < huge(0_i64)) then
      self%next_event_id = maximum_id + 1_i64
    else
      self%next_event_id = huge(0_i64)
    end if
    if (self%next_event_id < 1_i64 .or. self%next_event_id <= maximum_id) then
      error stop 'Restored outer-event next ID must exceed all active IDs.'
    end if
  end subroutine restore_outer_event_queue

  pure real(dp) function outer_event_signed_charge(self, species_id) result(charge)
    class(outer_event_queue_type), intent(in) :: self
    integer(i32), intent(in), optional :: species_id
    integer(i32) :: idx

    charge = 0.0_dp
    do idx = 1, self%count
      if (present(species_id)) then
        if (self%records(idx)%species_id /= species_id) cycle
      end if
      charge = charge + self%records(idx)%q*self%records(idx)%w
    end do
  end function outer_event_signed_charge

  pure real(dp) function outer_event_species_particle_number(self, species_id) result(particle_number)
    class(outer_event_queue_type), intent(in) :: self
    integer(i32), intent(in) :: species_id
    integer(i32) :: idx

    particle_number = 0.0_dp
    do idx = 1, self%count
      if (self%records(idx)%species_id == species_id) particle_number = particle_number + self%records(idx)%w
    end do
  end function outer_event_species_particle_number

  pure real(dp) function outer_event_photoelectron_number(self, photoelectron_species_id) result(particle_number)
    class(outer_event_queue_type), intent(in) :: self
    integer(i32), intent(in) :: photoelectron_species_id

    particle_number = self%species_particle_number(photoelectron_species_id)
  end function outer_event_photoelectron_number

  !> Canonical rank-local queue fingerprint components.
  !! The contract includes queue ordering, rank ownership, the next event ID, and
  !! every persisted field of every active event.  Reals use the checkpoint's
  !! round-trip-safe representation so a write/load cycle preserves the hash.
  subroutine outer_event_queue_fingerprint_components(self, mpi_rank, components)
    class(outer_event_queue_type), intent(in) :: self
    integer(i32), intent(in) :: mpi_rank
    integer(i64), intent(out) :: components(2)
    type(fingerprint_state_type) :: hash
    integer(i32) :: idx

    if (mpi_rank < 0_i32) error stop 'Outer-event queue fingerprint rank must be nonnegative.'
    call fingerprint_feed_string(hash, 'outer-event-queue-local-v1')
    call fingerprint_feed_i32(hash, mpi_rank)
    call fingerprint_feed_i64(hash, self%next_event_id)
    call fingerprint_feed_i32(hash, self%count)
    do idx = 1_i32, self%count
      call fingerprint_feed_event(hash, self%records(idx))
    end do
    components = [hash%a, hash%b]
  end subroutine outer_event_queue_fingerprint_components

  function outer_event_queue_fingerprint(self, mpi_rank) result(fingerprint)
    class(outer_event_queue_type), intent(in) :: self
    integer(i32), intent(in) :: mpi_rank
    character(len=16) :: fingerprint
    integer(i64) :: components(2)

    call self%fingerprint_components(mpi_rank, components)
    write (fingerprint, '(z8.8,z8.8)') components
  end function outer_event_queue_fingerprint

  !> Fold collectively summed local components into the restart-wide identity.
  function outer_event_queue_global_fingerprint(component_sums, world_size) result(fingerprint)
    integer(i64), intent(in) :: component_sums(2)
    integer(i32), intent(in) :: world_size
    character(len=16) :: fingerprint
    type(fingerprint_state_type) :: hash

    if (world_size <= 0_i32) error stop 'Outer-event queue fingerprint world size must be positive.'
    call fingerprint_feed_string(hash, 'outer-event-queue-global-v1')
    call fingerprint_feed_i32(hash, world_size)
    call fingerprint_feed_i64(hash, component_sums(1))
    call fingerprint_feed_i64(hash, component_sums(2))
    write (fingerprint, '(z8.8,z8.8)') hash%a, hash%b
  end function outer_event_queue_global_fingerprint

  pure logical function outer_event_queue_fingerprint_is_valid(fingerprint) result(valid)
    character(len=*), intent(in) :: fingerprint
    integer :: idx

    valid = len_trim(fingerprint) == 16
    if (.not. valid) return
    do idx = 1, 16
      select case (fingerprint(idx:idx))
      case ('0':'9', 'A':'F')
      case default
        valid = .false.
        return
      end select
    end do
  end function outer_event_queue_fingerprint_is_valid

  subroutine fingerprint_feed_event(hash, event)
    type(fingerprint_state_type), intent(inout) :: hash
    type(outer_event_record_type), intent(in) :: event

    call fingerprint_feed_i64(hash, event%event_id)
    call fingerprint_feed_real(hash, event%queued_time)
    call fingerprint_feed_real(hash, event%due_time)
    call fingerprint_feed_i32(hash, event%outcome)
    call fingerprint_feed_i32(hash, event%species_id)
    call fingerprint_feed_i32(hash, event%origin_rank)
    call fingerprint_feed_i32(hash, event%origin_batch)
    call fingerprint_feed_i64(hash, event%origin_particle_id)
    call fingerprint_feed_i32(hash, event%interface_face_index)
    call fingerprint_feed_real(hash, event%q)
    call fingerprint_feed_real(hash, event%m)
    call fingerprint_feed_real(hash, event%w)
    call fingerprint_feed_real_vector(hash, event%x)
    call fingerprint_feed_real_vector(hash, event%v)
  end subroutine fingerprint_feed_event

  subroutine fingerprint_feed_string(hash, value)
    type(fingerprint_state_type), intent(inout) :: hash
    character(len=*), intent(in) :: value
    integer :: idx

    call fingerprint_feed_byte(hash, len_trim(value))
    do idx = 1, len_trim(value)
      call fingerprint_feed_byte(hash, iachar(value(idx:idx)))
    end do
  end subroutine fingerprint_feed_string

  subroutine fingerprint_feed_i32(hash, value)
    type(fingerprint_state_type), intent(inout) :: hash
    integer(i32), intent(in) :: value
    character(len=32) :: encoded

    write (encoded, '(i0)') value
    call fingerprint_feed_string(hash, trim(encoded))
  end subroutine fingerprint_feed_i32

  subroutine fingerprint_feed_i64(hash, value)
    type(fingerprint_state_type), intent(inout) :: hash
    integer(i64), intent(in) :: value
    character(len=32) :: encoded

    write (encoded, '(i0)') value
    call fingerprint_feed_string(hash, trim(encoded))
  end subroutine fingerprint_feed_i64

  subroutine fingerprint_feed_real(hash, value)
    type(fingerprint_state_type), intent(inout) :: hash
    real(dp), intent(in) :: value
    character(len=32) :: encoded

    write (encoded, '(es26.17e3)') value
    call fingerprint_feed_string(hash, trim(adjustl(encoded)))
  end subroutine fingerprint_feed_real

  subroutine fingerprint_feed_real_vector(hash, values)
    type(fingerprint_state_type), intent(inout) :: hash
    real(dp), intent(in) :: values(:)
    integer :: idx

    call fingerprint_feed_i32(hash, int(size(values), i32))
    do idx = 1, size(values)
      call fingerprint_feed_real(hash, values(idx))
    end do
  end subroutine fingerprint_feed_real_vector

  subroutine fingerprint_feed_byte(hash, value)
    type(fingerprint_state_type), intent(inout) :: hash
    integer, intent(in) :: value

    hash%a = modulo(hash%a*fingerprint_multiplier_a + int(value, i64) + 1_i64, fingerprint_modulus)
    hash%b = modulo(hash%b*fingerprint_multiplier_b + int(value, i64) + 1_i64, fingerprint_modulus)
  end subroutine fingerprint_feed_byte

  subroutine ensure_outer_event_capacity(self, required)
    class(outer_event_queue_type), intent(inout) :: self
    integer(i32), intent(in) :: required
    type(outer_event_record_type), allocatable :: grown(:)
    integer(i32) :: old_capacity, new_capacity

    old_capacity = int(size(self%records), i32)
    if (required <= old_capacity) return
    if (old_capacity > huge(0_i32)/2_i32) then
      new_capacity = required
    else
      new_capacity = max(required, max(1_i32, 2_i32*old_capacity))
    end if
    allocate (grown(new_capacity))
    if (self%count > 0_i32) grown(1:self%count) = self%records(1:self%count)
    call move_alloc(grown, self%records)
  end subroutine ensure_outer_event_capacity

  subroutine validate_event_payload(event)
    type(outer_event_record_type), intent(in) :: event

    if (.not. ieee_is_finite(event%queued_time) .or. .not. ieee_is_finite(event%due_time)) then
      error stop 'Outer-event times must be finite.'
    end if
    if (event%queued_time < 0.0_dp .or. event%due_time <= event%queued_time) then
      error stop 'Outer-event flight interval must be finite and positive.'
    end if
    if (event%outcome /= outer_event_outcome_return .and. event%outcome /= outer_event_outcome_escape) then
      error stop 'Outer-event outcome must be return or escape.'
    end if
    if (event%species_id < 1_i32) error stop 'Outer-event species ID must be positive.'
    if (event%origin_rank < 0_i32 .or. event%origin_batch < 1_i32 .or. event%origin_particle_id < 1_i64) then
      error stop 'Outer-event origin metadata is invalid.'
    end if
    if (.not. ieee_is_finite(event%q)) error stop 'Outer-event charge must be finite.'
    if (.not. ieee_is_finite(event%m) .or. event%m <= 0.0_dp) then
      error stop 'Outer-event mass must be finite and positive.'
    end if
    if (.not. ieee_is_finite(event%w) .or. event%w <= 0.0_dp) then
      error stop 'Outer-event macro weight must be finite and positive.'
    end if
    if (.not. all(ieee_is_finite(event%x)) .or. .not. all(ieee_is_finite(event%v))) then
      error stop 'Outer-event phase-space coordinates must be finite.'
    end if
  end subroutine validate_event_payload

  pure logical function event_precedes(lhs, rhs) result(precedes)
    type(outer_event_record_type), intent(in) :: lhs, rhs

    precedes = lhs%due_time < rhs%due_time .or. &
               (lhs%due_time == rhs%due_time .and. lhs%event_id < rhs%event_id)
  end function event_precedes

  subroutine sort_events(events, count)
    type(outer_event_record_type), intent(inout) :: events(:)
    integer(i32), intent(in) :: count
    type(outer_event_record_type), allocatable :: work(:)
    integer(i32) :: width, left, middle, right, src_left, src_right, dst

    if (count < 2_i32) return
    allocate (work(count))
    width = 1_i32
    do while (width < count)
      left = 1_i32
      do while (left <= count)
        middle = min(left + width, count + 1_i32)
        right = min(left + 2_i32*width - 1_i32, count)
        src_left = left
        src_right = middle
        do dst = left, right
          if (src_left >= middle) then
            work(dst) = events(src_right)
            src_right = src_right + 1_i32
          else if (src_right > right) then
            work(dst) = events(src_left)
            src_left = src_left + 1_i32
          else if (event_precedes(events(src_right), events(src_left))) then
            work(dst) = events(src_right)
            src_right = src_right + 1_i32
          else
            work(dst) = events(src_left)
            src_left = src_left + 1_i32
          end if
        end do
        left = right + 1_i32
      end do
      events(1:count) = work
      if (width > count/2_i32) exit
      width = 2_i32*width
    end do
  end subroutine sort_events

  subroutine assert_unique_event_ids(events, count)
    type(outer_event_record_type), intent(in) :: events(:)
    integer(i32), intent(in) :: count
    integer(i64), allocatable :: ids(:), work(:)
    integer(i32) :: width, left, middle, right, src_left, src_right, dst, idx

    if (count < 2_i32) return
    allocate (ids(count), work(count))
    do idx = 1, count
      ids(idx) = events(idx)%event_id
    end do
    width = 1_i32
    do while (width < count)
      left = 1_i32
      do while (left <= count)
        middle = min(left + width, count + 1_i32)
        right = min(left + 2_i32*width - 1_i32, count)
        src_left = left
        src_right = middle
        do dst = left, right
          if (src_left >= middle) then
            work(dst) = ids(src_right)
            src_right = src_right + 1_i32
          else if (src_right > right) then
            work(dst) = ids(src_left)
            src_left = src_left + 1_i32
          else if (ids(src_right) < ids(src_left)) then
            work(dst) = ids(src_right)
            src_right = src_right + 1_i32
          else
            work(dst) = ids(src_left)
            src_left = src_left + 1_i32
          end if
        end do
        left = right + 1_i32
      end do
      ids = work
      if (width > count/2_i32) exit
      width = 2_i32*width
    end do
    do idx = 2, count
      if (ids(idx) == ids(idx - 1_i32)) error stop 'Restored outer-event IDs must be unique.'
    end do
  end subroutine assert_unique_event_ids

end module bem_outer_event_queue
