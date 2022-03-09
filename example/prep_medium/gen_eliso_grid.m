% Generate a low velocity basin model 

clear all;
close all;

%------------------------------------------------------------------------------
%-- x and y sampling
%------------------------------------------------------------------------------

nghost = 3;
nx = 126;
ny = 106;
dx = 100;
dy = 100;
x0 = 0.0 - nghost * dx;
y0 = 0.0 - nghost * dy;
Lx = nx * dx;
Ly = ny * dy;

x1d = [0 : nx-1] * dx + x0;
y1d = [0 : ny-1] * dy + y0;

%------------------------------------------------------------------------------
%-- velocity structures
%------------------------------------------------------------------------------

num_of_intfce = 3;
num_of_point_per_lay(1) = 11;
num_of_point_per_lay(2) = 21;
num_of_point_per_lay(3) = 101;
dz = 100;

[x2d, y2d] = meshgrid(x1d,y1d);

%-- backround 1D structure

%- roughly depth of each interface
dep     = [0, 2e3, 5e3];
Vp      = [ 1500, 3000, 5000];
Vs       = [ 1000, 1500, 2500 ];
den      = [ 1000, 2000, 3000 ];

%
%-- set inteface topo first
%

lay_elev = zeros(ny, nx, num_of_intfce);

%-- 1st: free surface
lay_elev(:,:,1) = 0.0;
%-- 2nd: top
num_circle_x = 3; %- how many circle along x
num_circle_y = 3; %- how many circle along y
lay_elev(:,:,2) =    cos(num_circle_x * x2d / Lx * 2 * pi)  ...
                  .* cos(num_circle_y * y2d / Ly * 2 * pi) ...
                  .* 5 * dx ...
                  - dep(2);
%-- 3rd: flat
lay_elev(:,:,3) = -dep(3);

%
%-- set points of each interface
%

for n = 1 : num_of_intfce
  point_elev{n} = zeros(ny, nx, num_of_point_per_lay(n));
  point_Vp{n} = zeros(ny, nx, num_of_point_per_lay(n));
  point_Vs{n} = zeros(ny, nx, num_of_point_per_lay(n));
  point_den{n} = zeros(ny, nx, num_of_point_per_lay(n));
end

%-- 1st: free surface
ilay = 1;
for j = 1 : ny
for i = 1 : nx
  %-- top and bot of this column
  elev_top = lay_elev(j,i,ilay);
  elev_bot = lay_elev(j,i,ilay+1);

  for ipt = 1 : num_of_point_per_lay(ilay)
    point_Vp  {ilay}(j,i,ipt) =  Vp(ilay) + ipt * 10;
    point_Vs  {ilay}(j,i,ipt) =  Vs(ilay) + ipt * 20;
    point_den {ilay}(j,i,ipt) =  den(ilay);
    point_elev{ilay}(j,i,ipt) = elev_top ...
          + (elev_bot-elev_top) / (num_of_point_per_lay(ilay)-1) *(ipt-1);
  end
end
end

%-- 2nd: 
ilay = 2;
for j = 1 : ny
for i = 1 : nx
  %-- top and bot of this column
  elev_top = lay_elev(j,i,ilay);
  elev_bot = lay_elev(j,i,ilay+1);

  for ipt = 1 : num_of_point_per_lay(ilay)
    point_Vp  {ilay}(j,i,ipt) =  Vp(ilay) - ipt * 10;
    point_Vs  {ilay}(j,i,ipt) =  Vs(ilay) - ipt * 10;
    point_den {ilay}(j,i,ipt) =  den(ilay);
    point_elev{ilay}(j,i,ipt) = elev_top ...
          + (elev_bot-elev_top) / (num_of_point_per_lay(ilay)-1) *(ipt-1);
  end
end
end

%-- 3rd: last layer
ilay = 3;
for j = 1 : ny
for i = 1 : nx
  %-- top and bot of this column
  elev_top = lay_elev(j,i,ilay);
  elev_bot = elev_top - (num_of_point_per_lay(ilay)-1) * dz;

  for ipt = 1 : num_of_point_per_lay(ilay)
    point_Vp  {ilay}(j,i,ipt) =  Vp(ilay);
    point_Vs  {ilay}(j,i,ipt) =  Vs(ilay);
    point_den {ilay}(j,i,ipt) =  den(ilay);
    point_elev{ilay}(j,i,ipt) = elev_top ...
          + (elev_bot-elev_top) / (num_of_point_per_lay(ilay)-1) *(ipt-1);
  end
end
end

%------------------------------------------------------------------------------
%-- plot
%------------------------------------------------------------------------------

figure;

for n = 1 : num_of_intfce
    surf( x2d, y2d,  ...
          point_elev{n}(:,:,1), ...
          point_Vp  {n}(:,:,1));
    hold on;

    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    %shading flat;
    shading faceted;
    %title('structure grid');
    colorbar;
    title('Vp at top of each layers');
end

%------------------------------------------------------------------------------
%-- create md3grd structure
%------------------------------------------------------------------------------

md.media_type = 'elastic_isotropic'
%  one_component, 
%  acoustic_isotropic, 
%  elastic_isotropic, 
%  elastic_vti_prem, elastic_vti_thomsen, elastic_vti_cij,
%  elastic_tti_thomsen, elastic_tti_bond,
%  elastic_aniso_cij

md.num_of_intfce = num_of_intfce

md.num_of_point_per_lay = num_of_point_per_lay;

md.nx = nx;
md.ny = ny;
md.dx = dx;
md.dy = dy;
md.x0 = x0;
md.y0 = y0;

%-- elastic iso
%md.elev = point_elev;
%md.Vp = point_Vp;
%md.Vs = point_Vs;
%md.density = point_den;
%-- permute if order is not [x,y,layer]
for n = 1 : num_of_intfce
  md.elev{n}    = permute(point_elev{n},[2,1,3]);
  md.density{n} = permute(point_den {n},[2,1,3]);
  md.Vp  {n}    = permute(point_Vp  {n},[2,1,3]);
  md.Vs  {n}    = permute(point_Vs  {n},[2,1,3]);
end

%------------------------------------------------------------------------------
%-- export
%------------------------------------------------------------------------------

fnm_ou = 'topolay_el_iso.md3grd'

md3grd_export(fnm_ou, md);

