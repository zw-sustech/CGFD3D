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

basin_x0 = x1d(50);
basin_y0 = y1d(53);
basin_e0 = 0.0 ;
basin_r0 = 15 * dx;

%------------------------------------------------------------------------------
%-- velocity structures
%------------------------------------------------------------------------------

num_of_layer = 3;

%-- backround 1D structure

%- roughly depth of each interface
dep     = [-1e3, basin_r0, 3e3, 10e3];
%dep     = [0, basin_r0, 3e3, 10e3];

Vp      = [ 1500, 3000, 5000, 8000];
Vp_grad = [ 0.2 ,  0.0, 0.0,  0 ];
Vp_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

den      = [ 1000, 2000, 3000, 3500 ];
den_grad = [ 0.2 ,  0.0, 0.0,  0 ];
den_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

%-- construct 3D structure,
%--   grad and pow are general 1D enough

[x2d, y2d] = meshgrid(x1d,y1d);

lay_elev = zeros(ny, nx, num_of_layer+1);
lay_Vp   = zeros(ny, nx, num_of_layer+1);
lay_den  = zeros(ny, nx, num_of_layer+1);


%-- 1st: free surface
lay_Vp  (:,:,1) =   Vp(1);
lay_den (:,:,1) =  den(1);
lay_elev(:,:,1) = -dep(1);

%-- 2nd: basin
lay_Vp  (:,:,2) = Vp(2);
lay_den (:,:,2) =den(2);

lay_elev(:,:,2) = lay_elev(:,:,1);
for j = 1 : ny
for i = 1 : nx
    r2d = sqrt((x1d(i) - basin_x0)^2 + (y1d(j) - basin_y0)^2);
    %if in basin
    if r2d <= basin_r0
      %-- (x-x0)^2 + (y-y0)^2 + (e-e0)^2 = r^2
      lay_elev(j,i,2) = - (sqrt( basin_r0^2 - r2d^2 ) + basin_e0 );
    end
end
end

% velocity perturbation
num_circle_x = 4; %- how many circle along x
num_circle_y = 3; %- how many circle along y
lay_Vp(:,:,2) =  (   cos(num_circle_x * x2d / Lx * 2*pi)  ...
                  .* cos(num_circle_y * y2d/Ly*2*pi) ...
                  .* 0.2 ...
                  + 1.0 ) ...
                  .* lay_Vp(:,:,2);
lay_den(:,:,2) = (    cos(num_circle_x * x2d / Lx * 2*pi)  ...
                  .* cos(num_circle_y * y2d/Ly*2*pi) ...
                  .* 0.2 ...
                  + 1.0 ) ...
                  .* lay_den(:,:,2);

%-- 3rd: topo
lay_Vp  (:,:,3) =  Vp(3);
lay_den (:,:,3) = den(3);

num_circle_x = 3; %- how many circle along x
num_circle_y = 1; %- how many circle along y

lay_elev(:,:,3) =    cos(num_circle_x * x2d / Lx * 2*pi)  ...
                  .* cos(num_circle_y * y2d/Ly*2*pi) ...
                  .* 5 * dx ...
                  - dep(3);

%-- 4rd: topo
lay_Vp  (:,:,4) =  Vp(4);
lay_den (:,:,4) = den(4);
lay_elev(:,:,4) = -dep(4);

%------------------------------------------------------------------------------
%-- plot
%------------------------------------------------------------------------------

figure;
for ilay = 1 : num_of_layer+1
    mesh(x2d, y2d, lay_elev(:,:,ilay));
    hold on;
    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    title('The build model');
end

%%%% 
%figure;
%drawmodel(501,501,300,-2500,0,-3000,10,10,10,'rho.dat',[],0,[]);

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
media_type = 'acoustic_isotropic'

fnm_ou = 'basin_ac_iso.md3lay'

fid = fopen(fnm_ou,'w');

%-- 1st: type
fprintf(fid, '%s\n',media_type);

%-- 2nd: number of layer
fprintf(fid, '%d\n', num_of_layer + 1);

%-- 3rd
fprintf(fid, '%d %d %f %f %f %f\n', nx, ny, x0, y0, dx, dy);

%-- rest
  for ilay = 1 : num_of_layer+1
      for j = 1 : ny
          for i = 1 : nx
              % elevation
            	fprintf(fid, '%g', lay_elev(j,i,ilay));
              % rho
            	fprintf(fid, ' %g %g %g', lay_den(j,i,ilay), den_grad(ilay), den_pow(ilay));
              % Vp
            	fprintf(fid, ' %g %g %g', lay_Vp(j,i,ilay), Vp_grad(ilay), Vp_pow(ilay));
              % return
            	fprintf(fid, '\n');
          end
      end
  end
fclose(fid);

