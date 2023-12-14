%------------------------------------------------------------------------------
%-- Example of generating source file using discrete stf values
%------------------------------------------------------------------------------

clear all;
close all;

%-- event name
src.evtnm = 'moment';

%-- number of point sources
src.number_of_source = 1;

%-- stf given by name or discrete values
%   0 analytic stf or 1 discrete values
src.stf_is_discrete = 1;

%for analytical, 2nd value is time length
src.stf_time_length = 0.0;
%for discrete,  2nd value is dt and 3rd is nt
src.stf_dt = 0.0002;
src.stf_nt = 30001;

%-- cmp,  1(force), 2(momoment), 3(force+moment)
src.cmp_fi_mij = 2;

%-- mechanism type: 0 by mij, 1 by angle + mu + D + A
src.mechansim_type = 0;

%-- location by indx or axis: 0 grid index, 1 coords
src.loc_coord_type = 1;

%-- 3rd dim is depth or not
src.loc_3dim = 0;

%-- coords
src.x_coord = [0 ];
src.y_coord = [0 ];
src.z_coord = [0 ];

%-- start time of each point source
src.t_start   = [0];

%-- stf 
T=0.1;
%-- mech
mxx = 0.0;
myy = 0.0;
mzz = 0.0;
myz = 0.0;
mxz = 0.0;
mxy = 1e18;

% single-force
Fxx = 0;
Fyy = 0;
Fzz = 1e18;


fc=1.5;
t0=1;
for is = 1 : src.number_of_source
  for it = 1 : src.stf_nt
      %- time for stf cal
      t = (it-1) * src.stf_dt;

      src.Mxx(it, is) = mxx * stf_moment(t,T);
      src.Myy(it, is) = myy * stf_moment(t,T);
      src.Mzz(it, is) = mzz * stf_moment(t,T);
      src.Myz(it, is) = myz * stf_moment(t,T);
      src.Mxz(it, is) = mxz * stf_moment(t,T);
      src.Mxy(it, is) = mxy * stf_moment(t,T);
      

%       src.Mxx(it, is) = mxx * stf_ricker(t,fc,t0);
%       src.Myy(it, is) = myy * stf_ricker(t,fc,t0);
%       src.Mzz(it, is) = mzz * stf_ricker(t,fc,t0);
%       src.Myz(it, is) = myz * stf_ricker(t,fc,t0);
%       src.Mxz(it, is) = mxz * stf_ricker(t,fc,t0);
%       src.Mxy(it, is) = mxy * stf_ricker(t,fc,t0);

%         src.Fx(it, is) = Fxx * stf_moment(t,T);
%         src.Fy(it, is) = Fyy * stf_moment(t,T);
%         src.Fz(it, is) = Fzz * stf_moment(t,T);

  end
end

%------------------------------------------------------------------------------
%-- export
%------------------------------------------------------------------------------

fnm = 'moment.src';

src_export(fnm, src);

% function
function S=stf_moment(t,T)
S = (t).*exp(-(t)/T)/(T*T);

end


