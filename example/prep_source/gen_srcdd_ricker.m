%------------------------------------------------------------------------------
%-- Example of generating nc source file using discrete ricker wavelet values
%------------------------------------------------------------------------------

clear all;
close all;

%-- number of point sources
src.number_of_source = 3;

%-- stf time
src.stf_t0 = 0.0;
src.stf_dt = 0.02;
src.stf_nt = 501;

%-- location by indx or axis: 0 grid index, 1 coords
src.loc_coord_type = 1
%-- 3rd dim is depth or not
src.loc_3dim = 1;

%-- coords
src.x_coord = [10e3, 14e3, 20e3 ];
src.y_coord = [10e3, 14e3, 20e3 ];
src.z_coord = [550.0,  550.0, 550.0 ];

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
      t         = (it-1) * src.stf_dt;
      src.t(it) = t;

      src.Mxx(it,is) = mxx * stf_ricker(t,ricker_fc(is),ricker_t0(is));
      src.Myy(it,is) = myy * stf_ricker(t,ricker_fc(is),ricker_t0(is));
      src.Mzz(it,is) = mzz * stf_ricker(t,ricker_fc(is),ricker_t0(is));
      src.Myz(it,is) = myz * stf_ricker(t,ricker_fc(is),ricker_t0(is));
      src.Mxz(it,is) = mxz * stf_ricker(t,ricker_fc(is),ricker_t0(is));
      src.Mxy(it,is) = mxy * stf_ricker(t,ricker_fc(is),ricker_t0(is));
  end
end

%------------------------------------------------------------------------------
%-- create to nc file
%--  use matlab native netcdf support
%------------------------------------------------------------------------------

fnm = 'event_3moment_srcdd.nc';

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

exit
