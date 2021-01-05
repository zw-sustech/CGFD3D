clc;clear;
output_dir='export/home/wangyh/CGFD3D-elastic/project/output';
%output_dir='/home/zhangw/work/cgfd_mpi/output.6.free'
%output_dir='/home/zhangw/work/cgfd_mpi/output.3.1src.discount'

nprocx = 1;
nprocy = 3;

%figure

i_indx = 49;
vname = 'Vx'

figure;
snap_fname = {};

  for ipx = 2    %%all in px2 process ,so....
  for ipy = 1 : nprocy
    snap_fname{ipx-1, ipy} = [output_dir, '/', ...
        'slicex_i', num2str(i_indx), ...
        '_px', num2str(ipx-1), '_py' num2str(ipy-1), ...
        '.nc'];
    
  end
  end
  
  %%
  tmpx = [];
  var = [];

for ipy = 1 : nprocy      
        tmp = ncread(snap_fname{ipx-1 ,ipy}, vname );
        tmpx = vertcat(tmpx, tmp);
end
        

  
%%
var = permute(tmpx,[2 1 3]);;
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
   caxis([-6,6]*1e-10);
  %caxis([-1,1]*1e-12);

  colorbar
  drawnow;
%   pause(0.5);

%     figure;
%    plot(squeeze(var(1,50,:)),'-*')
%    figure;
%   plot(squeeze(var(26,50,:)),'-*')
%  % plot(squeeze(var(63,:,43)),'-*')
% 
%   figure
%   plot(squeeze(var(30,50,:)),'-*')
  %figure
  %plot(squeeze(var00(63,23,:)),'-*')
end

