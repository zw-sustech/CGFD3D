%------------------------------------------------------------------------------
%-- Example of generating source file using stf function name
%------------------------------------------------------------------------------

clear all;
close all;

%-- event name
src.evtnm = 'event_3force_ricker';

%-- number of point sources
src.number_of_source = 3;

%-- stf given by name or discrete values
%   0 analytic stf or 1 discrete values
src.stf_is_discrete = 0;

%for analytical, 2nd value is time length
src.stf_time_length = 4.0;
%for discrete,  2nd value is dt and 3rd is nt
src.stf_dt = 0.0;
src.stf_nt = 1;

%-- cmp,  1(force), 2(momoment), 3(force+moment)
src.cmp_fi_mij = 0;

%-- mechanism type: 0 by mij, 1 by angle + mu + D + A
src.mechansim_type = 0;

%-- location by indx or axis: 0 grid index, 1 coords
src.loc_coord_type = 1

%-- 3rd dim is depth or not
src.loc_3dim = 0;

%-- coords
src.x_coord = [1000.0, 2000.0, 3000.0 ];
src.y_coord = [2200.0, 3050.0, 5000.0 ];
src.z_coord = [3200.0, 2300.0, 1000.0 ];

%-- stf
src.stf_coefs = zeros(10, src.number_of_source);

src.stf_name = ['ricker'; 'ricker'; 'ricker'];

ricker_fc = [2.0, 2.0, 2.0];
ricker_t0 = [0.5, 0.5, 0.5];
src.stf_coefs(1:2,:) = [ricker_fc; ricker_t0];

%-- start time of each point source
src.t_start   = [1.0,0.0,2.0];

% force_vector and moment_tensor
src.Fx = [1e16, 1e16, 1e16];
src.Fy = [1e16, 1e16, 1e16];
src.Fz = [1e15, 1e15, 1e15];

%- if cmp_fi_mij >= 1
%src.Mxx = [1e15, 1e15, 1e15];
%src.Myy = [1e15, 1e15, 1e15];
%src.Mzz = [1e15, 1e15, 1e15];
%src.Myz = [1e15, 1e15, 1e15];
%src.Mxz = [1e15, 1e15, 1e15];
%src.Mxy = [1e15, 1e15, 1e15];

%- if mechansim_type == 1
%src.strike = [80.0, 70.0, 50.0];
%src.dip    = [70.0, 60.0, 50.0];
%src.rake   = [30.0, 40.0, 50.0];
%src.mu     = [1e10, 1e10, 1e10]; % u is shear modulus
%src.D      = [1.0, 1.1, 1.2]; %unit is meter
%src.A     =  [1.0e8, 1.0e8, 1.0e8]; %unit is meter^2

%------------------------------------------------------------------------------
%-- export
%------------------------------------------------------------------------------

fnm = 'event_3force_ricker.src';

src_export(fnm, src);
