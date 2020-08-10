#ifndef SV_ELISO1ST_CURV_MACDRP_H
#define SV_ELISO1ST_CURV_MACDRP_H

void
sv_eliso1st_curv_macdrp_allstep(
    float *restrict w3d,  // wavefield
    float *restrict g3d,  // grid vars
    float *restrict m3d,  // medium vars
    // grid size
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    // boundary type
    int *restrict boundary_itype,
    // if abs
    int              abs_itype, //
    int    *restrict abs_num_of_layers, //
    size_t *restrict abs_indx, //
    size_t *restrict abs_coefs_facepos0, //
    float  *restrict abs_coefs, //
    size_t           abs_vars_size_per_level, //
    size_t *restrict abs_vars_volsiz, //
    size_t *restrict abs_vars_facepos0, //
    float  *restrict abs_vars,
    // if free surface
    float *matVx2Vz, //
    float *matVy2Vz, //
    // source term
    int num_of_force,
    int *restrict force_info, // num_of_force * 6 : si,sj,sk,start_pos_in_stf,start_it, end_it
    float *restrict force_vec_stf,
    int num_of_moment,
    int *restrict moment_info, // num_of_force * 6 : si,sj,sk,start_pos_in_rate,start_it, end_it
    float *restrict moment_ten_rate,
    // io
    int num_of_sta, int *restrict sta_loc_point, float *restrict sta_seismo,
    int num_of_snap, int *restrict snap_grid_indx, int *restrict snap_time_indx,
    // scheme
    int num_rk_stages, float *rk_a, float *rk_b, int num_of_pairs, 
    size_t fdx_max_len, size_t fdy_max_len, size_t fdz_max_len, size_t fdz_num_surf_lay,
    size_t ****pair_fdx_all_info, size_t ***pair_fdx_all_indx, float ***pair_fdx_all_coef,
    size_t ****pair_fdy_all_info, size_t ***pair_fdy_all_indx, float ***pair_fdy_all_coef,
    size_t ****pair_fdz_all_info, size_t ***pair_fdz_all_indx, float ***pair_fdz_all_coef,
    // time
    float dt, int nt_total, float t0,
    // mpi
    int myid, int *myid2, MPI_Comm comm,
    int qc_check_nan_num_of_step,
    const int verbose, // used for fprint qc
    char *out_dir);

void
sv_eliso1st_curv_macdrp_onestage(
    float *restrict w_cur, float *restrict rhs, 
    float *restrict g3d, float *restrict m3d,
    // grid size
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    int *restrict boundary_itype,
    int              abs_itype,
    int    *restrict abs_num_of_layers,
    size_t *restrict abs_indx,
    size_t *restrict abs_coefs_facepos0,
    float  *restrict abs_coefs,
    size_t           abs_vars_size_per_level,
    size_t *restrict abs_vars_volsiz,
    size_t *restrict abs_vars_facepos0, //
    float  *restrict abs_vars_cur,
    float  *restrict abs_vars_rhs,
    float *matVx2Vz, float *matVy2Vz, //
    // source term
    int num_of_force,
    int *restrict force_loc_point,
    float *restrict force_vec_value, // only for cur stage, size: num_of_force
    int num_of_moment,
    int *restrict moment_loc_point,
    float *restrict moment_ten_value,
    // include different order/stentil
    int fdx_max_len, int fdy_max_len, int fdz_max_len, int fdz_num_surf_lay,
    size_t **restrict fdx_all_info, size_t *restrict fdx_all_indx,float *restrict fdx_all_coef,
    size_t **restrict fdy_all_info, size_t *restrict fdy_all_indx,float *restrict fdy_all_coef,
    size_t **restrict fdz_all_info, size_t *restrict fdz_all_indx,float *restrict fdz_all_coef,
    const int myid, const int verbose);

void
sv_eliso1st_curv_macdrp_rhs_inner(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    const int myid, const int verbose);

void
sv_eliso1st_curv_macdrp_rhs_timg_z2(
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict jac3d, float *restrict slw3d,
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    const int myid, const int verbose);

void
sv_eliso1st_curv_macdrp_rhs_vlow_z2(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    float *restrict matVx2Vz, float *restrict matVy2Vz,
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_num_surf_lay, size_t fdz_max_len, size_t **restrict fdz_all_info,
    size_t *restrict fdz_all_indx, float *restrict fdz_all_coef,
    const int myid, const int verbose);

void
sv_eliso1st_curv_macdrp_rhs_cfspml(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict  mu3d, float *restrict slw3d,
    size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    int *restrict boundary_itype,
    int *restrict abs_num_of_layers,
    size_t *restrict abs_indx, // pml range of each face
    size_t *restrict abs_coefs_facepos0, // pos of 0 elem of each face
    float  *restrict abs_coefs,
    size_t *restrict abs_vars_volsiz, // size of single var in abs_blk_vars for each block
    size_t *restrict abs_vars_facepos0, // start pos of each blk in first level; other level similar
    float  *restrict abs_vars_cur, //
    float  *restrict abs_vars_rhs,
    const int myid, const int verbose);

void
sv_eliso1st_curv_macdrp_rhs_cfspml_timg_z2(
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict jac3d, float *restrict slw3d,
    size_t nk2, size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    int *restrict boundary_itype,
    int *restrict abs_num_of_layers,
    size_t *restrict abs_indx, // pml range of each face
    size_t *restrict abs_coefs_facepos0, // pos of 0 elem of each face
    float  *restrict abs_coefs,
    size_t *restrict abs_vars_volsiz, // size of single var in abs_blk_vars for each block
    size_t *restrict abs_vars_facepos0, // start pos of each blk in first level; other level similar
    float  *restrict abs_vars_cur, //
    float  *restrict abs_vars_rhs,
    const int myid, const int verbose);

void
sv_eliso1st_curv_macdrp_rhs_cfspml_vfree_z2(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    float *restrict matVx2Vz, float *restrict matVy2Vz,
    size_t nk2, size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    int *restrict boundary_itype,
    int *restrict abs_num_of_layers,
    size_t *restrict abs_indx, // pml range of each face
    size_t *restrict abs_coefs_facepos0, // pos of 0 elem of each face
    float  *restrict abs_coefs,
    size_t *restrict abs_vars_volsiz, // size of single var in abs_blk_vars for each block
    size_t *restrict abs_vars_facepos0, // start pos of each blk in first level; other level similar
    float  *restrict abs_vars_cur, //
    float  *restrict abs_vars_rhs,
    const int myid, const int verbose);

void
sv_eliso1st_curv_macdrp_rhs_src(
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict jac3d, float *restrict slw3d,
    size_t siz_line, size_t siz_slice,
    int num_of_force,
    int *restrict force_loc_point,
    float *restrict force_vec_value,
    int num_of_moment,
    int *restrict moment_loc_point,
    float *restrict moment_ten_value, // size: num_of_moment * 6
    const int myid, const int verbose);

void
sv_eliso1st_curv_macdrp_vel_dxy2dz(
    float *restrict g3d,
    float *restrict m3d,
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    float **restrict p_matVx2Vz,
    float **restrict p_matVy2Vz,
    const int myid, const int verbose);


#endif
