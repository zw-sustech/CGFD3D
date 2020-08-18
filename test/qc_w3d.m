
output_dir='/home/zhangw/work/cgfd_init_run/output';

px = 0;
py = 0;

for it = 0:10:200
  snap_fname = [output_dir, '/', ...
      'w3d_', num2str(px), '_' num2str(py), ...
      '_it', num2str(it), '.nc'];

  var = nc_varget(snap_fname,'Vz');

  pcolor(squeeze(var(32,:,:)));
  title(['it=',num2str(it)]);
  %%caxis([-1,1]*1e-25);
  colorbar
  drawnow;
end
