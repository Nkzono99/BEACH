#ifndef BEACH_FIELD_KERNEL_H
#define BEACH_FIELD_KERNEL_H

#ifdef __cplusplus
extern "C" {
#endif

#define BEACH_FIELD_KERNEL_ABI_MAJOR 1
#define BEACH_FIELD_KERNEL_ABI_MINOR 0

typedef void *beach_kernel_handle;

enum beach_kernel_status {
  BEACH_KERNEL_OK = 0,
  BEACH_KERNEL_INVALID_HANDLE = 1,
  BEACH_KERNEL_INVALID_ARGUMENT = 2,
  BEACH_KERNEL_NOT_READY = 3
};

enum beach_kernel_far_correction {
  BEACH_KERNEL_FAR_AUTO = 0,
  BEACH_KERNEL_FAR_NONE = 1,
  BEACH_KERNEL_FAR_RESERVED_2 = 2,
  BEACH_KERNEL_FAR_CACHED_KNEQ0 = 3
};

/*
 * Coordinate and vector arrays use xyz-interleaved storage:
 * values[3 * point_index + component]. Input arrays remain owned by the caller.
 * A handle returned by beach_kernel_create must be released exactly once with
 * beach_kernel_destroy.
 */
int beach_kernel_get_abi_version(int *major, int *minor);
int beach_kernel_get_build_info(char *buffer, int buffer_capacity,
                                int *text_length);
int beach_kernel_create(beach_kernel_handle *handle);
int beach_kernel_destroy(beach_kernel_handle handle);

int beach_kernel_set_periodic_cache(beach_kernel_handle handle,
                                    const char *path, int path_length,
                                    double tolerance);
int beach_kernel_get_periodic_cache_info(
    beach_kernel_handle handle, int *cache_hit, int *operator_build_count,
    char *fingerprint, int fingerprint_capacity, int *fingerprint_length,
    char *path, int path_capacity, int *path_length);

int beach_kernel_build(
    beach_kernel_handle handle, int source_count, const double *source_xyz,
    double theta, int leaf_max, int order, double softening,
    int use_periodic2, const int *periodic_axes_1based,
    const double *periodic_lengths, int image_layers, int far_correction,
    double ewald_alpha, int ewald_layers, const double *box_min,
    const double *box_max);
int beach_kernel_build_panel(
    beach_kernel_handle handle, int source_count, const double *vertex0_xyz,
    const double *vertex1_xyz, const double *vertex2_xyz, double theta,
    int leaf_max, int order, double softening, int use_periodic2,
    const int *periodic_axes_1based, const double *periodic_lengths,
    int image_layers, int far_correction, double ewald_alpha,
    int ewald_layers, const double *box_min, const double *box_max);
int beach_kernel_update_charges(beach_kernel_handle handle, int source_count,
                                const double *source_charges);

int beach_kernel_eval_e(beach_kernel_handle handle, int target_count,
                        const double *target_xyz, double *field_xyz);
int beach_kernel_eval_phi(beach_kernel_handle handle, int target_count,
                          const double *target_xyz, double *potential);
int beach_kernel_eval_e_direct(beach_kernel_handle handle, int target_count,
                               const double *target_xyz, double *field_xyz);
int beach_kernel_eval_phi_direct(beach_kernel_handle handle, int target_count,
                                 const double *target_xyz,
                                 double *potential);
int beach_kernel_force_on_charges(
    beach_kernel_handle handle, int target_count, const double *target_xyz,
    const double *target_charges, const double origin_xyz[3],
    double force_xyz[3], double torque_xyz[3]);

#ifdef __cplusplus
}
#endif

#endif
