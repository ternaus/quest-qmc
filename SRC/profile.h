#include <time.h>

/* perfomance counters */

#ifdef DQMC_PROFILE

#define PROFILE_MAX 10

#define PROFILE_ENABLE() profile_enabled = 1;

#define PROFILE_DISABLE() profile_enabled = 0;

#define PROFILE_BEGIN() if (profile_enabled) { \
  profile_count++; \
  profile_time[profile_count] = get_time(); \
  profile_flops[profile_count] = 0; }

#define PROFILE_END(i, n) if (profile_enabled) { \
  profile[i].count++; \
  profile[i].time += get_time() - profile_time[profile_count]; \
  profile[i].flops += n + profile_flops[profile_count]; \
  profile_flops[profile_count - 1] += n + profile_flops[profile_count]; \
  profile_count--; }

#define PROFILE_LWORK_BEGIN() PROFILE_BEGIN()
#define PROFILE_LWORK_END(i, n) PROFILE_END(i, n) 

enum {
  profile_dcopy, profile_daxpy, profile_dscal, profile_ddot, profile_dgemv, 
  profile_dgemm, profile_dtrmm, profile_dtrsm,
  profile_dsyev, profile_dgesv, profile_dgeqp3, profile_dorgqr, 
  profile_dormqr, profile_dgetrf, profile_dgetri,
  profile_dgetrs, profile_dgejsv, profile_dgerfsx, profile_dgeqrf, profile_dlarft, 
  profile_scalerow, profile_scalerowcol, 
  profile_normcol, profile_sort, profile_permute, profile_scalerowperm, profile_scalerowadd,
  profile_computeg, profile_computem, profile_swapg, 
  profile_transfer, profile_last
};

struct profile_item { 
  const char* name; int count; double time; double flops; 
};

extern profile_item profile[profile_last];

extern int profile_enabled, profile_count;
extern double profile_time[PROFILE_MAX];
extern double profile_flops[PROFILE_MAX];

extern "C" {
  void profile_enable_();
  void profile_disable_();
  void profile_begin_();
  void profile_end_(int *i, int *n);
  void profile_print_();
  void get_time_(double *t);
}

inline double get_time() {
  struct timespec tv;
  clock_gettime(CLOCK_MONOTONIC, &tv);
  return tv.tv_sec + tv.tv_nsec * 1e-9;
}

#else

#define PROFILE_ENABLE()
#define PROFILE_DISABLE()
#define PROFILE_BEGIN()
#define PROFILE_END(i, n)
#define PROFILE_LWORK_BEGIN()
#define PROFILE_LWORK_END(i, n) 

#endif
