#include <cstdio>

#include "profile.h"

#ifdef DQMC_PROFILE

profile_item profile[profile_last] = { 
  { "DCOPY    ", 0, 0, 0 }, 
  { "DAXPY    ", 0, 0, 0 }, 
  { "DSCAL    ", 0, 0, 0 }, 
  { "DDOT     ", 0, 0, 0 }, 
  { "DGEMV    ", 0, 0, 0 }, 
  { "DGEMM    ", 0, 0, 0 },
  { "DTRMM    ", 0, 0, 0 },
  { "DTRSM    ", 0, 0, 0 },
  { "DSYEV    ", 0, 0, 0 },
  { "DGESV    ", 0, 0, 0 },
  { "DGEQP3   ", 0, 0, 0 },
  { "DORGQR   ", 0, 0, 0 },
  { "DORMQR   ", 0, 0, 0 },
  { "DGETRF   ", 0, 0, 0 },
  { "DGETRI   ", 0, 0, 0 },
  { "DGETRS   ", 0, 0, 0 },
  { "DGEJSV   ", 0, 0, 0 },
  { "DGERFSX  ", 0, 0, 0 },
  { "DGEQRF   ", 0, 0, 0 },
  { "DLARFT   ", 0, 0, 0 },
  { "ScaleRow ", 0, 0, 0 },
  { "SRCol    ", 0, 0, 0 },
  { "NormCol  ", 0, 0, 0 },
  { "Sort     ", 0, 0, 0 },
  { "Permute  ", 0, 0, 0 },
  { "SRPerm   ", 0, 0, 0 },
  { "SRAdd    ", 0, 0, 0 },
  { "ComputeG ", 0, 0, 0 },
  { "ComputeM ", 0, 0, 0 },
  { "SwapG    ", 0, 0, 0 },
  { "Transfer ", 0, 0, 0 }
};

int profile_enabled = 0, profile_count = 0;
double profile_time[PROFILE_MAX];
double profile_flops[PROFILE_MAX];

void profile_enable_() { PROFILE_ENABLE(); }

void profile_disable_() { PROFILE_DISABLE(); }

void profile_begin_() { PROFILE_BEGIN(); }

void profile_end_(int *i, int *n) { PROFILE_END(*i, *n); }

void profile_print_()
{
  fflush(NULL);
  for (int i = 0; i < profile_last; i++) {
    if (profile[i].count) 
      printf(" %s %18i  %22.16f %10.2f\n", profile[i].name, profile[i].count, profile[i].time,
        profile[i].time > 0.0 ? profile[i].flops * 1e-6 / profile[i].time : 0.0);
  } 
}

void get_time_(double *t) { *t = get_time(); }

#endif
