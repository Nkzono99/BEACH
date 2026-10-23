!> 粒子注入APIを責務別internal moduleから再公開する互換facade。
module bem_injection
  use bem_injection_random, only: seed_rng, sample_uniform_positions, &
                                  sample_shifted_maxwell_velocities, init_random_beam_particles
  use bem_injection_geometry, only: compute_face_area_from_bounds
  use bem_injection_flux, only: compute_inflow_flux_from_drifting_maxwellian, &
                                compute_macro_particles_for_batch, compute_macro_particles_from_flux, &
                                sample_reservoir_face_particles, sample_reservoir_face_particle_counter
  use bem_photoelectron_injection, only: sample_photo_raycast_particles, cast_photo_ray, sample_photo_ray_slot, &
                                         photo_ray_macro_weight
  implicit none
  private

  public :: seed_rng
  public :: sample_uniform_positions
  public :: sample_shifted_maxwell_velocities
  public :: init_random_beam_particles
  public :: compute_inflow_flux_from_drifting_maxwellian
  public :: compute_face_area_from_bounds
  public :: compute_macro_particles_for_batch
  public :: compute_macro_particles_from_flux
  public :: sample_reservoir_face_particles
  public :: sample_reservoir_face_particle_counter
  public :: cast_photo_ray, sample_photo_ray_slot, photo_ray_macro_weight
  public :: sample_photo_raycast_particles

end module bem_injection
