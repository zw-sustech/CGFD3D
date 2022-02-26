% Draw seismic wavefield snapshot on multi cross-sections by using surf style
% Author:       Yuanhang Huo
% Email:        yhhuo@mail.ustc.edu.cn
% Affiliation:  University of Science and Technology of China
% Date:         2021.06.09

clear all;
addmypath
% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% which snapshot to plot
% profile 1
id{1}=1;
subs{1}=[20,1,1];      % start from index '1'
subc{1}=[1,-1,-1];     % '-1' to plot all points in this dimension
subt{1}=[1,1,1];
% profile 2
id{2}=1;
subs{2}=[1,1,30];      % start from index '1'
subc{2}=[-1,-1,1];     % '-1' to plot all points in this dimension
subt{2}=[1,1,1];
% profile 3
id{3}=1;
subs{3}=[1,30,1];      % start from index '1'
subc{3}=[-1,1,-1];     % '-1' to plot all points in this dimension
subt{3}=[1,1,1];

% variable and time to plot
varnm='Vz';
ns=1;
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



% locate snapshot and load coordinate
for i=1:length(subs)
    
    % locate snapshot
    snapinfo{i}=locate_snap(parfnm,id{i},'start',subs{i},'count',subc{i},'stride',subt{i},'snapdir',output_dir);
    % get coordinate data
    [x{i},y{i},z{i}]=gather_coord(snapinfo{i},'coorddir',output_dir);
    nx{i}=size(x{i},1);
    ny{i}=size(x{i},2);
    nz{i}=size(x{i},3);
    % coordinate unit
    str_unit='m';
    if flag_km
        x{i}=x{i}/1e3;
        y{i}=y{i}/1e3;
        z{i}=z{i}/1e3;
        str_unit='km';
    end

end

% figure plot
hid=figure;
set(hid,'BackingStore','on');

% snapshot show
for nlayer=ns:nt:ne
    
    for i=1:length(subs)
        
        [v{i},t]=gather_snap(snapinfo{i},nlayer,varnm,'snapdir',output_dir);
        
        % show time
        if i==1
            disp([ '  draw ' num2str(nlayer) 'th time step (t=' num2str(t) ')']);
        end
        
        % show snapshot
        if flag_emlast
            sid{i}=surf(squeeze(permute(x{i},[2 1 3])), ...
                squeeze(permute(y{i},[2 1 3])), ...
                squeeze(permute(z{i},[2 1 3])), ...
                squeeze(permute(v{i},[2 1 3])));
        else
            sid{i}=surf(flipdim(squeeze(permute(x{i},[2 1 3])),3), ...
                flipdim(squeeze(permute(y{i},[2 1 3])),3), ...
                flipdim(squeeze(permute(z{i},[2 1 3])),3), ...
                flipdim(squeeze(permute(v{i},[2 1 3])),3));
        end
        
        hold on;
        
    end
    
    hold off;
        
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


