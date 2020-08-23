output_dir='/home/zhangw/work/cgfd_mpi/output'

nprocx = 2;
nprocy = 1;

%figure

%for it = 0:10:1500
for it = 110 %0:10:1500
  snap_fname = [output_dir, '/', ...
      'w3d_0_0', ...
      '_it', num2str(it), '.nc'];
  var00 = nc_varget(snap_fname,'Vz');

  snap_fname = [output_dir, '/', ...
      'w3d_1_0', ...
      '_it', num2str(it), '.nc'];
  var10 = nc_varget(snap_fname,'Vz');

  %var = [var00(1:66,1:106,1:53);var10(1:66,1:106,4:56)];
  var = cat(3,var00(1:66,1:106,1:53),var10(1:66,1:106,4:56));

  %pcolor(squeeze(var(:,:,43)));
  pcolor(squeeze(var(:,23,:)));
  %pcolor(squeeze(var(59,:,:)));

  %pcolor(squeeze(var00(59,:,:)));
  %pcolor(squeeze(var10(59,:,:)));

  %title(['it=',num2str(it)]);
  title(['it=',num2str(it), ' of ', output_dir]);
  caxis([-1,1]*1e-10);
  %caxis([-1,1]*1e-12);
  colorbar
  drawnow;
end

