clc;
clear all;
% close all;

addmypath
% -------------------------- parameters input -------------------------- %
% file and path name
parfnm='../project/test.json';
output_dir='../project/output';
PG_dir = output_dir;
% get free surface x y coords
%subs is start index, subc is counts, subt is step.
subs=[1,1,1];      % start from index '1'
subc=[-1,-1,1];     % '-1' to plot all points in this dimension
subt=[1,1,1];

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
scl_daspect =[1 1 1];
clrmp       = 'parula';

% variable to plot
% 'PGV', 'PGVh', 'PGVx', 'PGVy', 'PGVz', 'PGA', 'PGAh', 'PGAx',
% 'PGAy','PGAz', 'PGD', 'PGDh', 'PGDx', 'PGDy','PGDz'
varnm='PGA';

% ---------------------------------------------------------------------- %

% load grid coordinate
coordinfo=locate_coord(parfnm,'start',subs,'count',subc,'stride',subt,'coorddir',output_dir);
[x,y,z]=gather_coord(coordinfo,'coorddir',output_dir);
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

% readl all PG_V_A_D nc files
PG_prefix='PG_V_A_D';
PG_list=dir([PG_dir,'/',PG_prefix,'*.nc']);

for i=1:length(PG_list)    
    PG_nm=[PG_dir,'/',PG_list(i).name];
    topoid = nc_attget(PG_nm,nc_global,'coords_of_mpi_topo');
    px(i) = topoid(1);
    py(i) = topoid(2);
    counts = nc_attget(PG_nm,nc_global,'count_index_of_physical_points');
    ni = counts(1);
    nj = counts(2);
    gstart = nc_attget(PG_nm,nc_global,'global_index_of_first_physical_points');
    gni1 = gstart(1)+1;
    gnj1 = gstart(2)+1;
    subs1 = [3,3]; 
    subc1 = [nj,ni]; 
    subt1 = [1,1]; 
    fnm_PG = [PG_dir,'/',PG_prefix,'_px',num2str(px(i)),'_py',num2str(py(i)),'.nc']
    V(gnj1:gnj1+nj-1,gni1:gni1+ni-1)=nc_varget(fnm_PG,varnm,subs1,subc1,subt1);
end
% transpose to (x,y)
V = V';
%%

% figure plot
hid=figure;
set(hid,'BackingStore','on');
pcolor(x,y,V);
xlabel(['X axis (',str_unit,')']);
ylabel(['Y axis (',str_unit,')']);
shading flat;
set(gcf,'color','white','renderer','painters');
colormap( 'jet' );
colorbar();
% caxis([0,3]);
title(varnm);

% save and print figure
if flag_print==1
    width= 500;
    height=500;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    fnm_out=[varnm];
    print(gcf,[fnm_out '.png'],'-dpng');
end
