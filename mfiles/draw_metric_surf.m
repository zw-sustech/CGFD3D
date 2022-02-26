% Draw a cross-section of metric by using surf style
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

% which metric profile to plot
subs=[10,30,1];     % start from index '1'
subc=[-1,-1,1];     % '-1' to plot all points in this dimension
subt=[2,1,2];

% variable to plot
% 'jac', 'xi_x', 'xi_y', 'xi_z', 'eta_x', 'eta_y', 'eta_z',
% 'zeta_x', 'zeta_y', 'zeta_z'
varnm='zeta_x';

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_clb    = 1;
flag_title  = 1;
scl_daspect = [1 1 1];
clrmp       = 'parula';
% ---------------------------------------------------------------------- %



% locate metric data
metricinfo=locate_metric(parfnm,'start',subs,'count',subc,'stride',subt,'metricdir',output_dir);
% get coordinate data
[x,y,z]=gather_coord(metricinfo,'coorddir',output_dir);
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

% gather metric data
v=gather_metric(metricinfo,varnm,'metricdir',output_dir);

% figure plot
hid=figure;
set(hid,'BackingStore','on');

% metric show
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


