%------------------------------------------------------------------------------
%-- Example of generating source file using discrete stf values
%------------------------------------------------------------------------------

clear all;
close all;

%-- event name
src.evtnm = 'event_3moment_discrete';

%-- number of point sources
src.number_of_source = 3;

%-- stf given by name or discrete values
%   0 analytic stf or 1 discrete values
src.stf_is_discrete = 1;

%for analytical, 2nd value is time length
src.stf_time_length = 0.0;
%for discrete,  2nd value is dt and 3rd is nt
src.stf_dt = 0.02;
src.stf_nt = 51;

%-- cmp,  1(force), 2(momoment), 3(force+moment)
src.cmp_fi_mij = 2;

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

%-- start time of each point source
src.t_start   = [1.0,0.0,2.0];

%-- stf 
ricker_fc = [2.0, 2.0, 2.0];
ricker_t0 = [0.5, 0.5, 0.5];
%-- mech
mxx = 1e16;
myy = 1e16;
mzz = 1e16;
myz = 0.0;
mxz = 0.0;
mxy = 0.0;

for is = 1 : src.number_of_source
  for it = 1 : src.stf_nt
      %- time for stf cal
      t = (it-1) * src.stf_dt;

      src.Mxx(it, is) = mxx * stf_ricker(t,ricker_fc(is),ricker_t0(is));
      src.Myy(it, is) = myy * stf_ricker(t,ricker_fc(is),ricker_t0(is));
      src.Mzz(it, is) = mzz * stf_ricker(t,ricker_fc(is),ricker_t0(is));
      src.Myz(it, is) = myz * stf_ricker(t,ricker_fc(is),ricker_t0(is));
      src.Mxz(it, is) = mxz * stf_ricker(t,ricker_fc(is),ricker_t0(is));
      src.Mxy(it, is) = mxy * stf_ricker(t,ricker_fc(is),ricker_t0(is));
  end
end

%------------------------------------------------------------------------------
%-- export
%------------------------------------------------------------------------------

fnm = 'event_3moment_discrete.src';

src_export(fnm, src);

