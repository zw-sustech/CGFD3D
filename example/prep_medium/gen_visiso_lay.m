% Generate a low velocity basin model 

clear all;
close all;

%------------------------------------------------------------------------------
%-- x and y sampling
%------------------------------------------------------------------------------

nghost = 3;
nx = 126;
ny = 106;
dx = 300;
dy = 300;
x0 = 0.0 - nghost * dx;
y0 = 0.0 - nghost * dy;
Lx = nx * dx;
Ly = ny * dy;

x1d = [0 : nx-1] * dx + x0;
y1d = [0 : ny-1] * dy + y0;


%------------------------------------------------------------------------------
%-- velocity structures
%------------------------------------------------------------------------------

num_of_layer = 3;

%-- backround 1D structure

%- roughly depth of each interface
dep     = [0, 1.5e3, 3e3, 10e3];

Vp      = [ 1700, 3000, 5000, 8000];
Vp_grad = [ 0.2 ,  0.0, 0.0,  0 ];
Vp_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

Vs       = [ 1000, 1500, 2500, 3000 ];
Vs_grad  = [ 0.2 ,  0.0, 0.0,  0 ];
Vs_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

den      = [ 1000, 2000, 3000, 3500 ];
den_grad = [ 0.2 ,  0.0, 0.0,  0 ];
den_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

Qp      = [ 60, 70, 90 ,150];
Qp_grad = [ 0.2 ,  0.0, 0.0,  0 ];
Qp_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

Qs      = [ 30, 45, 50 ,80];
Qs_grad = [ 0.2 ,  0.0, 0.0,  0 ];
Qs_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

%-- construct 3D structure,
%--   grad and pow are general 1D enough

[x2d, y2d] = meshgrid(x1d,y1d);

%-- calculate elevation of each interface
lay_elev = zeros(ny, nx, num_of_layer+1);

lay_Vp      = zeros(ny, nx, num_of_layer+1);
lay_Vp_grad = zeros(ny, nx, num_of_layer+1);
lay_Vp_pow  = zeros(ny, nx, num_of_layer+1);

lay_Vs      = zeros(ny, nx, num_of_layer+1);
lay_Vs_grad = zeros(ny, nx, num_of_layer+1);
lay_Vs_pow  = zeros(ny, nx, num_of_layer+1);

lay_den  = zeros(ny, nx, num_of_layer+1);
lay_den_grad = zeros(ny, nx, num_of_layer+1);
lay_den_pow  = zeros(ny, nx, num_of_layer+1);

lay_Qp      = zeros(ny, nx, num_of_layer+1);
lay_Qp_grad = zeros(ny, nx, num_of_layer+1);
lay_Qp_pow  = zeros(ny, nx, num_of_layer+1);

lay_Qs      = zeros(ny, nx, num_of_layer+1);
lay_Qs_grad = zeros(ny, nx, num_of_layer+1);
lay_Qs_pow  = zeros(ny, nx, num_of_layer+1);

%-- 1st: free surface
lay_elev(:,:,1) = -dep(1);

lay_Vp     (:,:,1) = Vp     (1);
lay_Vp_grad(:,:,1) = Vp_grad(1);
lay_Vp_pow (:,:,1) = Vp_pow (1);

lay_Vs     (:,:,1) = Vs     (1);
lay_Vs_grad(:,:,1) = Vs_grad(1);
lay_Vs_pow (:,:,1) = Vs_pow (1);

lay_den     (:,:,1) =       den(1);
lay_den_grad(:,:,1) =  den_grad(1);
lay_den_pow (:,:,1) =  den_pow (1);

lay_Qp     (:,:,1) = Qp     (1);
lay_Qp_grad(:,:,1) = Qp_grad(1);
lay_Qp_pow (:,:,1) = Qp_pow (1);

lay_Qs     (:,:,1) = Qs     (1);
lay_Qs_grad(:,:,1) = Qs_grad(1);
lay_Qs_pow (:,:,1) = Qs_pow (1);

%-- 2nd: basin
lay_Vp     (:,:,2) = Vp     (2);
lay_Vp_grad(:,:,2) = Vp_grad(2);
lay_Vp_pow (:,:,2) = Vp_pow (2);

lay_Vs     (:,:,2) = Vs     (2);
lay_Vs_grad(:,:,2) = Vs_grad(2);
lay_Vs_pow (:,:,2) = Vs_pow (2);

lay_den     (:,:,2) =       den(2);
lay_den_grad(:,:,2) =  den_grad(2);
lay_den_pow (:,:,2) =  den_pow (2);

lay_Qp     (:,:,2) = Qp     (2);
lay_Qp_grad(:,:,2) = Qp_grad(2);
lay_Qp_pow (:,:,2) = Qp_pow (2);

lay_Qs     (:,:,2) = Qs     (2);
lay_Qs_grad(:,:,2) = Qs_grad(2);
lay_Qs_pow (:,:,2) = Qs_pow (2);

lay_elev(:,:,2) = -dep(2);

%-- 3rd: topo
lay_Vp     (:,:,3) = Vp     (3);
lay_Vp_grad(:,:,3) = Vp_grad(3);
lay_Vp_pow (:,:,3) = Vp_pow (3);

lay_Vs     (:,:,3) = Vs     (3);
lay_Vs_grad(:,:,3) = Vs_grad(3);
lay_Vs_pow (:,:,3) = Vs_pow (3);

lay_den     (:,:,3) =       den(3);
lay_den_grad(:,:,3) =  den_grad(3);
lay_den_pow (:,:,3) =  den_pow (3);

lay_Qp     (:,:,3) = Qp     (3);
lay_Qp_grad(:,:,3) = Qp_grad(3);
lay_Qp_pow (:,:,3) = Qp_pow (3);

lay_Qs     (:,:,3) = Qs     (3);
lay_Qs_grad(:,:,3) = Qs_grad(3);
lay_Qs_pow (:,:,3) = Qs_pow (3);

lay_elev(:,:,3) = - dep(3);

%-- 4rd: topo
lay_Vp     (:,:,4) = Vp     (4);
lay_Vp_grad(:,:,4) = Vp_grad(4);
lay_Vp_pow (:,:,4) = Vp_pow (4);

lay_Vs     (:,:,4) = Vs     (4);
lay_Vs_grad(:,:,4) = Vs_grad(4);
lay_Vs_pow (:,:,4) = Vs_pow (4);

lay_den     (:,:,4) =       den(4);
lay_den_grad(:,:,4) =  den_grad(4);
lay_den_pow (:,:,4) =  den_pow (4);

lay_Qp     (:,:,4) = Qp     (4);
lay_Qp_grad(:,:,4) = Qp_grad(4);
lay_Qp_pow (:,:,4) = Qp_pow (4);

lay_Qs     (:,:,4) = Qs     (4);
lay_Qs_grad(:,:,4) = Qs_grad(4);
lay_Qs_pow (:,:,4) = Qs_pow (4);

lay_elev(:,:,4) = -dep(4);

%------------------------------------------------------------------------------
%-- plot
%------------------------------------------------------------------------------
%%
figure;
for ilay = 1 : num_of_layer+1
    surf(x2d, y2d,lay_elev(:,:,ilay), lay_Qp(:,:,ilay));
    hold on;
    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    title('The build model');
    shading interp
    colorbar
end
%%
%%%% 
%figure;
%drawmodel(501,501,300,-2500,0,-3000,10,10,10,'rho.dat',[],0,[]);

%------------------------------------------------------------------------------
%-- create md3lay structure
%------------------------------------------------------------------------------

md.media_type = 'elastic_isotropic'
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

  md.Vs     {n} = permute(lay_Vs     (:,:,n),[2,1]);
  md.Vs_coef{n} = permute(lay_Vs_grad(:,:,n),[2,1]);
  md.Vs_pow {n} = permute(lay_Vs_pow (:,:,n),[2,1]);
end

%------------------------------------------------------------------------------
%-- export
%------------------------------------------------------------------------------

fnm_ou = 'topolay_el_iso.md3lay'

md3lay_export(fnm_ou, md);

%------------------------------------------------------------------------------
%-- create md3lay structure
%------------------------------------------------------------------------------

md.media_type = 'one_component'
%  one_component, 
%  acoustic_isotropic, 
%  elastic_isotropic, 
%  elastic_vti_prem, elastic_vti_thomsen, elastic_vti_cij,
%  elastic_tti_thomsen, elastic_tti_bond,
%  elastic_aniso_cij

%-- elastic iso
%md.elev = point_elev;
%md.Vp = point_Vp;
%md.Vs = point_Vs;
%md.density = point_den;
%-- permute if order is not [x,y,layer]
for n = 1 : md.num_of_intfce
  md.val     {n} = permute(lay_Qp     (:,:,n),[2,1]);
  md.val_coef{n} = permute(lay_Qp_grad(:,:,n),[2,1]);
  md.val_pow {n} = permute(lay_Qp_pow (:,:,n),[2,1]);
end

%------------------------------------------------------------------------------
%-- export
%------------------------------------------------------------------------------

fnm_ou = 'topolay_vis_iso_Qp.md3lay'

md3lay_export(fnm_ou, md);
%------------------------------------------------------------------------------
%-- create md3lay structure
%------------------------------------------------------------------------------

md.media_type = 'one_component'
%  one_component, 
%  acoustic_isotropic, 
%  elastic_isotropic, 
%  elastic_vti_prem, elastic_vti_thomsen, elastic_vti_cij,
%  elastic_tti_thomsen, elastic_tti_bond,
%  elastic_aniso_cij

%-- elastic iso
%md.elev = point_elev;
%md.Vp = point_Vp;
%md.Vs = point_Vs;
%md.density = point_den;
%-- permute if order is not [x,y,layer]
for n = 1 : md.num_of_intfce
  md.val     {n} = permute(lay_Qs     (:,:,n),[2,1]);
  md.val_coef{n} = permute(lay_Qs_grad(:,:,n),[2,1]);
  md.val_pow {n} = permute(lay_Qs_pow (:,:,n),[2,1]);
end

%------------------------------------------------------------------------------
%-- export
%------------------------------------------------------------------------------

fnm_ou = 'topolay_vis_iso_Qs.md3lay'

md3lay_export(fnm_ou, md);