
%output_dir='/home/zhangw/work/cgfd_init_run/output';
output_dir='/home/zhangw/work/cgfd_arc/05pml/output'

px = 0;
py = 0;

for it = 1:10:150
  snap_fname = [output_dir, '/', ...
      'w3d_', num2str(px), '_' num2str(py), ...
      '_it', num2str(it), '.nc'];

  var = nc_varget(snap_fname,'Vz');

  %pcolor(squeeze(var(51,:,:)));
  pcolor(squeeze(var(:,:,41)));
  %pcolor(squeeze(var(:,41,:)));
  title(['it=',num2str(it)]);
  %%caxis([-1,1]*1e-25);
  colorbar
  drawnow;
end
