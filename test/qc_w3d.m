
%output_dir='/home/zhangw/work/cgfd_ac/00/output';
%output_dir='/home/zhangw/work/cgfd_stg_ac/00/output';
output_dir='/home/zhangw/work/cgfd_stg_ac/01/output';

px = 0; py = 0;
%px = 1; py = 0;
%px = 1; py = 1;

%for it = 50 %1:10:500
%for it = 63 %- px=1 begin
for it = 10:10:500
  snap_fname = [output_dir, '/', ...
      'w3d_', num2str(px), '_' num2str(py), ...
      '_it', num2str(it), '.nc'];

  %var = nc_varget(snap_fname,'Vx');
  %var = nc_varget(snap_fname,'Vy');
  %var = nc_varget(snap_fname,'Vz');
  %var = nc_varget(snap_fname,'Txx');
  %var = nc_varget(snap_fname,'Tyy');
  %var = nc_varget(snap_fname,'Tzz');
  %var = nc_varget(snap_fname,'Tyz');
  %var = nc_varget(snap_fname,'Txz');
  %var = nc_varget(snap_fname,'Txy');

  var = nc_varget(snap_fname,'P');

  %pcolor(squeeze(var(51+2,:,:)));
  pcolor(squeeze(var(:,:,41+2)));
  %pcolor(squeeze(var(:,41,:)));
  title(['it=',num2str(it)]);
  %caxis([-1,1]*1e-25);
  %caxis([-1,1]*1e-35);
  colorbar
  drawnow;
end
