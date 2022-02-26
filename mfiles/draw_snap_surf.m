% Draw seismic wavefield snapshot by using surf style
% Author:       Yuanhang Huo
% Email:        yhhuo@mail.ustc.edu.cn
% Affiliation:  University of Science and Technology of China
% Date:         2021.05.31

clear all;
addmypath
% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% which snapshot to plot
id=1;
subs=[1,1,50];      % start from index '1'
subc=[-1,-1,1];     % '-1' to plot all points in this dimension
subt=[1,1,1];

% variable and time to plot
varnm='Vx';
ns=50;
ne=500;
nt=50;

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_light  = 0;
savegif = 1;
% scl_caxis=[-1.0 1.0];
filename1 = ['Vx.gif'];
scl_daspect =[1 1 1];
clrmp       = 'jetwr';
taut=0.5;
% ---------------------------------------------------------------------- %



% load snapshot data
snapinfo=locate_snap(parfnm,id,'start',subs,'count',subc,'stride',subt,'snapdir',output_dir);
% get coordinate data
[x,y,z]=gather_coord(snapinfo,'coorddir',output_dir);
nx=size(x,1);
ny=size(x,2);
nz=size(x,3);
% coordinate unit
str_unit='m';
if flag_km
   x=x/1e3;
   y=y/1e3;
   z=z/1e3;
   str_unit='km';
end

% figure plot
hid=figure;
set(hid,'BackingStore','on');

% snapshot show
for nlayer=ns:nt:ne
    
    [v,t]=gather_snap(snapinfo,nlayer,varnm,'snapdir',output_dir);
    
    disp([ '  draw ' num2str(nlayer) 'th time step (t=' num2str(t) ')']);
    
    if flag_emlast
        sid=surf(squeeze(permute(x,[2 1 3])), ...
            squeeze(permute(y,[2 1 3])), ...
            squeeze(permute(z,[2 1 3])), ...
            squeeze(permute(v,[2 1 3])));
    else
        sid=surf(flipdim(squeeze(permute(x,[2 1 3])),3), ...
            flipdim(squeeze(permute(y,[2 1 3])),3), ...
            flipdim(squeeze(permute(z,[2 1 3])),3), ...
            flipdim(squeeze(permute(v,[2 1 3])),3));
    end
    
    xlabel(['X axis (' str_unit ')']);
    ylabel(['Y axis (' str_unit ')']);
    zlabel(['Z axis (' str_unit ')']);

    set(gca,'layer','top');
    set(gcf,'color','white','renderer','painters');

    % axis image
    % shading interp;
    shading flat;
    % colorbar range/scale
    if exist('scl_caxis')
        caxis(scl_caxis);
    end
    % axis daspect
    if exist('scl_daspect')
        daspect(scl_daspect);
    end
    % colormap and colorbar
    if exist('clrmp')
        colormap(clrmp);
    end
    colorbar('vert');
    if flag_light
        view(-40,35);
        set(gca,'box','off');
        camlight(0,10,'local');
        lighting phong;
    end
    
    % title
    titlestr=['Snapshot of ' varnm ' at ' ...
              '{\fontsize{12}{\bf ' ...
              num2str((t),'%7.3f') ...
              '}}s'];
    title(titlestr);
    
    drawnow;
    pause(taut);
    %save gif
    if savegif
      im=frame2im(getframe(gcf));
      [imind,map]=rgb2ind(im,256);
      if nlayer==ns
        imwrite(imind,map,filename1,'gif','LoopCount',Inf,'DelayTime',0.5);
      else
        imwrite(imind,map,filename1,'gif','WriteMode','append','DelayTime',0.5);
      end
    end
    % save and print figure
    if flag_print==1
        width= 500;
        height=500;
        set(gcf,'paperpositionmode','manual');
        set(gcf,'paperunits','points');
        set(gcf,'papersize',[width,height]);
        set(gcf,'paperposition',[0,0,width,height]);
        fnm_out=[varnm '_ndim_',num2str(nlayer,'%5.5i')];
        print(gcf,[fnm_out '.png'],'-dpng');
    end
    
end


