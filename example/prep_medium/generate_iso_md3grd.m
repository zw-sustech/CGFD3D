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
num_circle_y = 1; %- how many circle along y
lay_elev(:,:,2) =    cos(num_circle_x * x2d / Lx * 2*pi)  ...
                  .* cos(num_circle_y * y2d/Ly*2*pi) ...
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

%==============================================================================
%-- write file
%==============================================================================

% first line of 3lay header is media_type, can take:  
%  one_component, 
%  acoustic_isotropic, 
%  elastic_isotropic, 
%  elastic_vti_prem, elastic_vti_thomsen, elastic_vti_cij,
%  elastic_tti_thomsen, elastic_tti_bond,
%  elastic_aniso_cij
media_type = 'elastic_isotropic'

fnm_ou = 'topolay_el_iso.md3grd'

fid = fopen(fnm_ou,'w');

%-- 1st: type
fprintf(fid, '%s\n',media_type);

%-- 2nd: number of layer
fprintf(fid, '%d\n', num_of_intfce);

%-- 3rd: number of points of each layer
for n = 1 : num_of_intfce
   	fprintf(fid, ' %d', num_of_point_per_lay(n)); 
end
fprintf(fid, '\n'); 

%-- 4rd
fprintf(fid, '%d %d %f %f %f %f\n', nx, ny, x0, y0, dx, dy);

%-- rest
  for ilay = 1 : num_of_intfce
      for k = 1 : num_of_point_per_lay(ilay)
      for j = 1 : ny
          for i = 1 : nx
              % elevation
            	fprintf(fid, '%g',  point_elev{ilay}(j,i,k));
              % rho
            	fprintf(fid, ' %g', point_den{ilay}(j,i,k));
              % Vp
            	fprintf(fid, ' %g', point_Vp{ilay}(j,i,k));
              % Vs
            	fprintf(fid, ' %g', point_Vs{ilay}(j,i,k));
              % return
            	fprintf(fid, '\n');
          end
      end
      end
  end
fclose(fid);

