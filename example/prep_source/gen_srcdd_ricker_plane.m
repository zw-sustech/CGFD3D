%------------------------------------------------------------------------------
%-- Example of generating nc source file using discrete ricker wavelet values
%------------------------------------------------------------------------------

clear all;
close all;

nx = 300
ny = 300
x0 = 0
y0 = 0
dx = 100
dy = 100

%-- number of point sources
src.number_of_source = nx * ny;

%-- stf time
src.stf_t0 = 0.0;
src.stf_dt = 0.02;
src.stf_nt = 501;

%-- location by indx or axis: 0 grid index, 1 coords
src.loc_coord_type = 1
%-- 3rd dim is depth or not
src.loc_3dim = 1;

%-- coords
src.x_coord = zeros(1,src.number_of_source);
src.y_coord = zeros(1,src.number_of_source);
src.z_coord = zeros(1,src.number_of_source);
src.Mxx     = zeros(src.stf_nt, src.number_of_source);
src.Myy     = zeros(src.stf_nt, src.number_of_source);
src.Mzz     = zeros(src.stf_nt, src.number_of_source);
src.Myz     = zeros(src.stf_nt, src.number_of_source);
src.Mxz     = zeros(src.stf_nt, src.number_of_source);
src.Mxy     = zeros(src.stf_nt, src.number_of_source);

src.stf_fc  = zeros(1,src.number_of_source);
src.stf_t0  = zeros(1,src.number_of_source);

for j = 1 : ny
for i = 1 : nx
  n = i + (j-1) * nx;
  src.x_coord(n) = x0 + (i-1) * dx;
  src.y_coord(n) = y0 + (j-1) * dy;
  src.z_coord(n) = 3e3;
  src.stf_fc(n) = 2.0;
  src.stf_t0(n) = 0.5;
end
end

%-- mech
mxx = 1e16;
myy = 1e16;
mzz = 1e16;
myz = 0.0;
mxz = 0.0;
mxy = 0.0;

for it = 1 : src.stf_nt
    %- time for stf cal
    t         = (it-1) * src.stf_dt;
    src.t(it) = t;
end

for is = 1 : src.number_of_source
    is
    src.Mxx(:,is) = mxx * stf_ricker(src.t,src.stf_fc(is),src.stf_t0(is));
    src.Myy(:,is) = myy * stf_ricker(src.t,src.stf_fc(is),src.stf_t0(is));
    src.Mzz(:,is) = mzz * stf_ricker(src.t,src.stf_fc(is),src.stf_t0(is));
    src.Myz(:,is) = myz * stf_ricker(src.t,src.stf_fc(is),src.stf_t0(is));
    src.Mxz(:,is) = mxz * stf_ricker(src.t,src.stf_fc(is),src.stf_t0(is));
    src.Mxy(:,is) = mxy * stf_ricker(src.t,src.stf_fc(is),src.stf_t0(is));
end

%------------------------------------------------------------------------------
%-- create to nc file
%--  use matlab native netcdf support
%------------------------------------------------------------------------------

fnm = 'event_plane_srcdd.nc';

disp(['writing ',fnm]);

%-- CLOBBER: overwrite exitsing file
ncid = netcdf.create(fnm,'CLOBBER');

%-- define time and n dim
tdimid = netcdf.defDim(ncid,'time',src.stf_nt);
ndimid = netcdf.defDim(ncid,'number',src.number_of_source);

%-- define time vars
tid = netcdf.defVar(ncid,'time','NC_FLOAT',[tdimid]);

%-- define coord vars
xid = netcdf.defVar(ncid,'x','NC_FLOAT',[ndimid]);
yid = netcdf.defVar(ncid,'y','NC_FLOAT',[ndimid]);
zid = netcdf.defVar(ncid,'z','NC_FLOAT',[ndimid]);

%-- define vars
mxxid = netcdf.defVar(ncid,'Mxx_rate','NC_FLOAT',[tdimid,ndimid]);
myyid = netcdf.defVar(ncid,'Myy_rate','NC_FLOAT',[tdimid,ndimid]);
mzzid = netcdf.defVar(ncid,'Mzz_rate','NC_FLOAT',[tdimid,ndimid]);
myzid = netcdf.defVar(ncid,'Myz_rate','NC_FLOAT',[tdimid,ndimid]);
mxzid = netcdf.defVar(ncid,'Mxz_rate','NC_FLOAT',[tdimid,ndimid]);
mxyid = netcdf.defVar(ncid,'Mxy_rate','NC_FLOAT',[tdimid,ndimid]);

%-- put att
%netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'), ...
%              'location_is_axis', src.loc_coord_type,'NC_INT');
%netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'), ...
%              'z_is_depth', src.loc_3dim,'NC_INT');
netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'), ...
              'location_is_axis', int32(src.loc_coord_type));
netcdf.putAtt(ncid,netcdf.getConstant('GLOBAL'), ...
              'z_is_depth', int32(src.loc_3dim));

%-- end def
netcdf.endDef(ncid);

%-- write time var
netcdf.putVar(ncid,tid,src.t);

%-- write coord var
netcdf.putVar(ncid,xid,src.x_coord);
netcdf.putVar(ncid,yid,src.y_coord);
netcdf.putVar(ncid,zid,src.z_coord);

%-- write mij var
netcdf.putVar(ncid,mxxid,src.Mxx);
netcdf.putVar(ncid,myyid,src.Myy);
netcdf.putVar(ncid,mzzid,src.Mzz);
netcdf.putVar(ncid,myzid,src.Myz);
netcdf.putVar(ncid,mxzid,src.Mxz);
netcdf.putVar(ncid,mxyid,src.Mxy);

%-- close
netcdf.close(ncid);

%exit
