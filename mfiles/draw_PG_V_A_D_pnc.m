clc;
clear all;
% close all;

addmypath
% -------------------------- parameters input -------------------------- %

nc_file = '~/work/cgfd3d/test6/output6/PG_V_A_D.nc'

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
scl_daspect =[1 1 1];
clrmp       = 'parula';

% variable to plot
% 'PGV', 'PGVh', 'PGVx', 'PGVy', 'PGVz', 'PGA', 'PGAh', 'PGAx',
% 'PGAy','PGAz', 'PGD', 'PGDh', 'PGDx', 'PGDy','PGDz'
varnm='PGVz';
%varnm='PGA';

% ---------------------------------------------------------------------- %

vinfo = ncinfo(nc_file, varnm);
dim_size = vinfo.Size

% load grid coordinate
x = ncread(nc_file, 'x');
y = ncread(nc_file, 'y');
z = ncread(nc_file, 'z');

% read var
V = ncread(nc_file, varnm);

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

%pcolor(x,y,V);
surf(x,y,z,V);

xlabel(['X axis (',str_unit,')']);
ylabel(['Y axis (',str_unit,')']);
zlabel(['Z axis (',str_unit,')']);

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
