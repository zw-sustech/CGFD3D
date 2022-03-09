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
lay_Vp      = zeros(ny, nx, num_of_layer+1);
lay_Vp_grad = zeros(ny, nx, num_of_layer+1);
lay_Vp_pow  = zeros(ny, nx, num_of_layer+1);
lay_den      = zeros(ny, nx, num_of_layer+1);
lay_den_grad = zeros(ny, nx, num_of_layer+1);
lay_den_pow  = zeros(ny, nx, num_of_layer+1);

%-- 1st: free surface
lay_elev(:,:,1) = -dep(1);

lay_Vp     (:,:,1) = Vp     (1);
lay_Vp_grad(:,:,1) = Vp_grad(1);
lay_Vp_pow (:,:,1) = Vp_pow (1);

lay_den     (:,:,1) =       den(1);
lay_den_grad(:,:,1) =  den_grad(1);
lay_den_pow (:,:,1) =  den_pow (1);

%-- 2nd: basin
lay_Vp     (:,:,2) = Vp     (2);
lay_Vp_grad(:,:,2) = Vp_grad(2);
lay_Vp_pow (:,:,2) = Vp_pow (2);

lay_den     (:,:,2) =       den(2);
lay_den_grad(:,:,2) =  den_grad(2);
lay_den_pow (:,:,2) =  den_pow (2);

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
lay_Vp     (:,:,3) = Vp     (3);
lay_Vp_grad(:,:,3) = Vp_grad(3);
lay_Vp_pow (:,:,3) = Vp_pow (3);

lay_den     (:,:,3) =       den(3);
lay_den_grad(:,:,3) =  den_grad(3);
lay_den_pow (:,:,3) =  den_pow (3);


num_circle_x = 3; %- how many circle along x
num_circle_y = 1; %- how many circle along y

lay_elev(:,:,3) =    cos(num_circle_x * x2d / Lx * 2*pi)  ...
                  .* cos(num_circle_y * y2d/Ly*2*pi) ...
                  .* 5 * dx ...
                  - dep(3);

%-- 4rd: topo
lay_elev(:,:,4) = -dep(4);

lay_Vp     (:,:,4) = Vp     (4);
lay_Vp_grad(:,:,4) = Vp_grad(4);
lay_Vp_pow (:,:,4) = Vp_pow (4);

lay_den     (:,:,4) =       den(4);
lay_den_grad(:,:,4) =  den_grad(4);
lay_den_pow (:,:,4) =  den_pow (4);

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

%------------------------------------------------------------------------------
%-- create md3lay structure
%------------------------------------------------------------------------------

md.media_type = 'acoustic_isotropic'
%  one_component, 
%  acoustic_isotropic, 
%  elastic_isotropic, 
%  elastic_vti_prem, elastic_vti_thomsen, elastic_vti_cij,
%  elastic_tti_thomsen, elastic_tti_bond,
%  elastic_aniso_cij

md.num_of_intfce = num_of_layer + 1;

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
for n = 1 : md.num_of_intfce
  md.elev{n}    = permute(lay_elev(:,:,n),[2,1]);

  md.density     {n} = permute(lay_den     (:,:,n),[2,1]);
  md.density_coef{n} = permute(lay_den_grad(:,:,n),[2,1]);
  md.density_pow {n} = permute(lay_den_pow (:,:,n),[2,1]);

  md.Vp     {n} = permute(lay_Vp     (:,:,n),[2,1]);
  md.Vp_coef{n} = permute(lay_Vp_grad(:,:,n),[2,1]);
  md.Vp_pow {n} = permute(lay_Vp_pow (:,:,n),[2,1]);
end

%------------------------------------------------------------------------------
%-- export
%------------------------------------------------------------------------------

fnm_ou = 'basin_ac_iso.md3lay'

md3lay_export(fnm_ou, md);

