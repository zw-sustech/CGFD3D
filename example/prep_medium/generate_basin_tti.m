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
%dep     = [-1e3, basin_r0, 3e3, 10e3];
dep     = [0, basin_r0, 3e3, 10e3];

Vp      = [ 1500, 3000, 5000, 8000];
Vp_grad = [ 0.2 ,  0.0, 0.0,  0 ];
Vp_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

Vs       = [ 1000, 1500, 2500, 3000 ];
Vs_grad  = [ 0.2 ,  0.0, 0.0,  0 ];
Vs_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

den      = [ 1000, 2000, 3000, 3500 ];
den_grad = [ 0.2 ,  0.0, 0.0,  0 ];
den_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

epsilon      = [ 0.5,   0.1, 0.2, 0.3];
epsilon_grad = [ 0.2 ,  0.0, 0.0,  0 ];
epsilon_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

gamma      = [ 0.4,   0.2, 0.1, 0.3];
gamma_grad = [ 0.2 ,  0.0, 0.0,  0 ];
gamma_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

delta      = [ 0.5,   0.1, 0.2, 0.3];
delta_grad = [ 0.2 ,  0.0, 0.0,  0 ];
delta_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

azi      = [ 45,   90, 0, 0];
azi_grad = [ 0.0 ,  0.0, 0.0,  0 ];
azi_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

dip      = [ 45,   90, 0.0, 0.0];
dip_grad = [ 0.0 ,  0.0, 0.0,  0 ];
dip_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

%-- construct 3D structure,
%--   grad and pow are general 1D enough

[x2d, y2d] = meshgrid(x1d,y1d);

%-- calculate elevation of each interface
lay_elev = zeros(ny, nx, num_of_layer+1);
lay_Vp   = zeros(ny, nx, num_of_layer+1);
lay_Vs   = zeros(ny, nx, num_of_layer+1);
lay_den  = zeros(ny, nx, num_of_layer+1);
lay_epsilon = zeros(ny, nx, num_of_layer+1);
lay_gamma   = zeros(ny, nx, num_of_layer+1);
lay_delta   = zeros(ny, nx, num_of_layer+1);
lay_azi   = zeros(ny, nx, num_of_layer+1);
lay_dip   = zeros(ny, nx, num_of_layer+1);

%-- 1st: free surface
lay_Vp     (:,:,1) =      Vp(1);
lay_Vs     (:,:,1) =      Vs(1);
lay_den    (:,:,1) =     den(1);
lay_epsilon(:,:,1) = epsilon(1);
lay_gamma  (:,:,1) = gamma  (1);
lay_delta  (:,:,1) = delta  (1);
lay_azi    (:,:,1) =   azi  (1);
lay_dip    (:,:,1) =   dip  (1);
lay_elev   (:,:,1) =    -dep(1);

%-- 2nd: basin
lay_Vp     (:,:,2) =      Vp(2);
lay_Vs     (:,:,2) =      Vs(2);
lay_den    (:,:,2) =     den(2);
lay_epsilon(:,:,2) = epsilon(2);
lay_gamma  (:,:,2) = gamma  (2);
lay_delta  (:,:,2) = delta  (2);
lay_azi    (:,:,2) =   azi  (2);
lay_dip    (:,:,2) =   dip  (2);
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

%-- 3rd: topo
lay_Vp     (:,:,3) =      Vp(3);
lay_Vs     (:,:,3) =      Vs(3);
lay_den    (:,:,3) =     den(3);
lay_epsilon(:,:,3) = epsilon(3);
lay_gamma  (:,:,3) = gamma  (3);
lay_delta  (:,:,3) = delta  (3);
lay_azi    (:,:,3) =   azi  (3);
lay_dip    (:,:,3) =   dip  (3);

num_circle_x = 3; %- how many circle along x
num_circle_y = 1; %- how many circle along y

lay_elev(:,:,3) =    cos(num_circle_x * x2d / Lx * 2*pi)  ...
                  .* cos(num_circle_y * y2d/Ly*2*pi) ...
                  .* 5 * dx ...
                  - dep(3);

%-- 4rd: topo
lay_Vp     (:,:,4) =      Vp(4);
lay_Vs     (:,:,4) =      Vs(4);
lay_den    (:,:,4) =     den(4);
lay_epsilon(:,:,4) = epsilon(4);
lay_gamma  (:,:,4) = gamma  (4);
lay_delta  (:,:,4) = delta  (4);
lay_azi    (:,:,4) =   azi  (4);
lay_dip    (:,:,4) =   dip  (4);
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
media_type = 'elastic_tti_thomsen'

fnm_ou = 'basin_el_tti.md3lay'

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
              % Vs
            	fprintf(fid, ' %g %g %g', lay_Vs(j,i,ilay), Vs_grad(ilay), Vs_pow(ilay));
              % epsilon
            	fprintf(fid, ' %g %g %g', lay_epsilon(j,i,ilay), epsilon_grad(ilay), epsilon_pow(ilay));
              % delta
            	fprintf(fid, ' %g %g %g', lay_delta(j,i,ilay), delta_grad(ilay), delta_pow(ilay));
              % gamma
            	fprintf(fid, ' %g %g %g', lay_gamma(j,i,ilay), gamma_grad(ilay), gamma_pow(ilay));
              % azimuth
            	fprintf(fid, ' %g %g %g', lay_azi(j,i,ilay), azi_grad(ilay), azi_pow(ilay));
              % dip
            	fprintf(fid, ' %g %g %g', lay_dip(j,i,ilay), dip_grad(ilay), dip_pow(ilay));
              % return
            	fprintf(fid, '\n');
          end
      end
  end
fclose(fid);

