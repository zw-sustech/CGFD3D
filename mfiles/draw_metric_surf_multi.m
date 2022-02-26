% Draw multi cross-sections of metric by using surf style
% Author:       Yuanhang Huo
% Email:        yhhuo@mail.ustc.edu.cn
% Affiliation:  University of Science and Technology of China
% Date:         2021.06.06

clear all;
addmypath
% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';

% metric profiles to plot
% profile 1
subs{1}=[50,1,1];      % start from index '1'
subc{1}=[1,-1,-1];     % '-1' to plot all points in this dimension
subt{1}=[1,1,1];
% profile 2
subs{2}=[1,1,30];      % start from index '1'
subc{2}=[-1,-1,1];     % '-1' to plot all points in this dimension
subt{2}=[1,1,1];
% profile 3
subs{3}=[1,50,1];      % start from index '1'
subc{3}=[-1,1,-1];     % '-1' to plot all points in this dimension
subt{3}=[1,1,1];

% variable to plot
% 'jac', 'xi_x', 'xi_y', 'xi_z', 'eta_x', 'eta_y', 'eta_z',
% 'zeta_x', 'zeta_y', 'zeta_z'
varnm='jac';

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_clb    = 1;
flag_title  = 1;
scl_daspect = [1 1 1];
clrmp       = 'parula';
% ---------------------------------------------------------------------- %



% figure plot
hid=figure;
set(hid,'BackingStore','on');

for i=1:length(subs)
    
    % locate metric
    metricinfo{i}=locate_metric(parfnm,'start',subs{i},'count',subc{i},'stride',subt{i},'metricdir',output_dir);
    % get coordinate data
    [x{i},y{i},z{i}]=gather_coord(metricinfo{i},'coorddir',output_dir);
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
    
    % gather metric data
    v{i}=gather_metric(metricinfo{i},varnm,'metricdir',output_dir);
    
    % metric show
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

xlabel(['X axis (' str_unit ')']);
ylabel(['Y axis (' str_unit ')']);
zlabel(['Z axis (' str_unit ')']);

set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');

% shading
% shading interp;
shading flat;
% colorbar range/scale
if exist('scl_caxis','var')
    caxis(scl_caxis);
end
% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis tight
% colormap and colorbar
if exist('clrmp')
    colormap(clrmp);
end
if flag_clb
    cid=colorbar;
end

% title
if flag_title
    title(varnm,'interpreter','none');
end

% save and print figure
if flag_print
    width= 500;
    height=500;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    print(gcf,[varnm '.png'],'-dpng');
end


