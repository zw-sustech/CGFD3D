% Draw seismic wavefield snapshot by using pcolor style
% Author:       Yuanhang Huo
% Email:        yhhuo@mail.ustc.edu.cn
% Affiliation:  University of Science and Technology of China
% Date:         2021.05.31

clear all;

% -------------------------- parameters input -------------------------- %
% file and path name
%parfnm='./project/test.json';
%output_dir='./project/output';
%parfnm='/home/zhangw/work/cgfd_cart/13nc/test.json'
%output_dir='/home/zhangw/work/cgfd_cart/13nc/output'
parfnm='/home/zhangw/work/cgfd_ac/00/test.json'
output_dir='/home/zhangw/work/cgfd_ac/00/output'

% which snapshot to plot
id=1;

%-- z slice
subs=[1,1,53];      % start from index '1'
subc=[-1,-1,1];     % '-1' to plot all points in this dimension
subt=[1,1,1];

%-- y slice
%subs=[1,41,1];      % start from index '1'
%subc=[-1,1,-1];     % '-1' to plot all points in this dimension
%subt=[1,1,1];

%-- x slice
%subs=[41,1,1];      % start from index '1'
%subc=[1,-1,-1];     % '-1' to plot all points in this dimension
%subt=[1,1,1];

% variable and time to plot
varnm='Vx';
ns=1;
ne=500;
nt=50;
%ns=2;
%ne=500;
%nt=2;

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
scl_daspect =[1 1 1];
clrmp       = 'parula';
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
    
    if nx==1
        if flag_emlast
            sid=pcolor((flipud(permute(squeeze(y),[2 1]))), ...
                (flipud(permute(squeeze(z),[2 1]))), ...
                (flipud(permute(squeeze(v),[2 1]))));
        else
            sid=pcolor((permute(squeeze(y),[2 1])), ...
                (permute(squeeze(z),[2 1])), ...
                (permute(squeeze(v),[2 1])));
        end
        xlabel(['Y axis (' str_unit ')']);
        ylabel(['Z axis (' str_unit ')']);
        
    elseif ny==1
        if flag_emlast
            sid=pcolor((flipud(permute(squeeze(x),[2 1]))), ...
                (flipud(permute(squeeze(z),[2 1]))), ...
                (flipud(permute(squeeze(v),[2 1]))));
        else
            sid=pcolor((permute(squeeze(x),[2 1])), ...
                (permute(squeeze(z),[2 1])), ...
                (permute(squeeze(v),[2 1])));
        end
        xlabel(['X axis (' str_unit ')']);
        ylabel(['Z axis (' str_unit ')']);
        
    else
        if flag_emlast
            sid=pcolor((flipud(permute(squeeze(x),[2 1]))), ...
                (flipud(permute(squeeze(y),[2 1]))), ...
                (flipud(permute(squeeze(v),[2 1]))));
        else
            sid=pcolor((permute(squeeze(x),[2 1])), ...
                (permute(squeeze(y),[2 1])), ...
                (permute(squeeze(v),[2 1])));
        end
        xlabel(['X axis (' str_unit ')']);
        ylabel(['Y axis (' str_unit ')']);
    end
    
    set(gca,'layer','top');
    set(gcf,'color','white','renderer','painters');

    % axis image
    % shading
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
    
    %title
    titlestr=['Snapshot of ' varnm ' at ' ...
              '{\fontsize{12}{\bf ' ...
              num2str((t),'%7.3f') ...
              '}}s'];
    title(titlestr);
    
    drawnow;
    pause(taut);
    
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


