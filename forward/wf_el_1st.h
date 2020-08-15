#ifndef WF_EL_1ST_H
#define WF_EL_1ST_H

#define WF_EL_1ST_NVAR    9

#define WF_EL_1ST_SEQ_Vx  0
#define WF_EL_1ST_SEQ_Vy  1
#define WF_EL_1ST_SEQ_Vz  2
#define WF_EL_1ST_SEQ_Txx 3
#define WF_EL_1ST_SEQ_Tyy 4
#define WF_EL_1ST_SEQ_Tzz 5
#define WF_EL_1ST_SEQ_Txz 6
#define WF_EL_1ST_SEQ_Tyz 7
#define WF_EL_1ST_SEQ_Txy 8

// cal index of var, not finished
#define M_WF_EL_1ST_IND(icmp,it,istage,nt,num_stage) \
  ((icmp) * (nt) * (num_stage) + (it) * (num_stage) + (istage))

void 
wf_el_1st_init_vars(
    size_t siz_volume,
    int number_of_levels,
    int *number_of_vars,
    float  **p_w3d,
    size_t **p_w3d_pos,
    char  ***p_w3d_name);

void
wf_el_1st_check_value(float *restrict w, size_t siz_volume);

static inline float *
wf_el_1st_getptr_Vx(float *restrict w, size_t siz_volume)
{
  return (w + WF_EL_1ST_SEQ_Vx * siz_volume);
}

static inline float *
wf_el_1st_getptr_Vz(float *restrict w, size_t siz_volume)
{
  return (w + WF_EL_1ST_SEQ_Vz * siz_volume);
}

#endif
