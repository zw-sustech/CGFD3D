output_dir='/home/zhangw/work/cgfd_io/output'
%output_dir='/home/zhangw/work/cgfd_mpi/output.6.free'
%output_dir='/home/zhangw/work/cgfd_mpi/output.3.1src.discount'

nprocx = 1;
nprocy = 1;

%figure

snap_name = 'volume_vel'
i_indx = 49;
vname = 'Vy'

figure;

  for ipx = 1 : nprocx
  for ipy = 1 : nprocy
    snap_fname = [output_dir, '/', ...
        snap_name, ...
        '_px', num2str(ipx-1), '_py' num2str(ipy-1), ...
        '.nc'];
  end
  end

  var = nc_varget(snap_fname, vname);

%for it = 1:10:1500
for it = 300 %0:10:1500
  %figure
  pcolor(squeeze(var(it,1,:,:)));

  %pcolor(squeeze(var(:,23,:)));
  %title(['it=',num2str(it)]);
  title(['it=',num2str(it), ' of ', output_dir]);
  %caxis([-1,1]*1e-8);
  caxis([-1,1]*1e-10);
  %caxis([-1,1]*1e-12);
  colorbar
  drawnow;
  pause(0.5);

  %figure;
  %%plot(squeeze(var(59,23,:)),'-*')
  %plot(squeeze(var(63,23,:)),'-*')
  %%plot(squeeze(var(63,:,43)),'-*')

  %figure
  %plot(squeeze(var10(63,23,:)),'-*')
  %figure
  %plot(squeeze(var00(63,23,:)),'-*')
end

