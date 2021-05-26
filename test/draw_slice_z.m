output_dir='export/home/wangyh/CGFD3D-elastic/project/output';
%output_dir='/home/zhangw/work/cgfd_mpi/output.6.free'
%output_dir='/home/zhangw/work/cgfd_mpi/output.3.1src.discount'

nprocx = 3;
nprocy = 3;

%figure

i_indx = 59;
vname = 'Vx'

figure;
snap_fname = {};

  for ipx = 1 : nprocx
  for ipy = 1 : nprocy
    snap_fname{ipx, ipy} = [output_dir, '/', ...
        'slicez_k', num2str(i_indx), ...
        '_px', num2str(ipx-1), '_py' num2str(ipy-1), ...
        '.nc'];
    
  end
  end
  
  %%
  tmpx = [];
  var = [];
  
for ipx = 1 : nprocx
        tmpy = [];
for ipy = 1 : nprocy      
        tmp = ncread(snap_fname{ipx ,ipy}, vname );
        tmpy = horzcat(tmpy, tmp);
end
        tmpx = vertcat(tmpx, tmpy);
end
  
%%
var = tmpx;
%   var = nc_varget(snap_fname, vname);
%   var = ncread(snap_fname, vname );

for it  = 10:10:500
   grid off;
    pcolor(squeeze(var(:,:,it)));
    shading interp;
  %pcolor(squeeze(var(:,23,:)));
  %title(['it=',num2str(it)]);
  title(['it=',num2str(it), ' of ', output_dir]);
  %caxis([-1,1]*1e-8);
   caxis([-1,1]*1e-9);
  %caxis([-1,1]*1e-12);
  
grid off;
 colorbar
  drawnow;
  grid off;
%   pause(0.5);

%   figure;
 % plot(squeeze(var(59,23,:)),'-*')
  %plot(squeeze(var(63,23,:)),'-*')
 % plot(squeeze(var(63,:,43)),'-*')

  %figure
%  plot(squeeze(var(52,52,:)),'-*')
  %figure
  %plot(squeeze(var00(63,23,:)),'-*')
end
