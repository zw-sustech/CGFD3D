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
den_grad = [0, 0, 0, 0];
den_pow  = [0, 0, 0, 0];
C_grad = [0, 0, 0, 0];
C_pow  = [0, 0, 0, 0];

%-- construct 3D structure,
%--   grad and pow are general 1D enough

[x2d, y2d] = meshgrid(x1d,y1d);

%-- calculate elevation of each interface
lay_elev = zeros(ny, nx, num_of_layer+1);
lay_den      = zeros(ny, nx, num_of_layer+1);
lay_den_grad = zeros(ny, nx, num_of_layer+1);
lay_den_pow  = zeros(ny, nx, num_of_layer+1);

lay_C    = zeros(21, ny, nx, num_of_layer+1);
lay_C_grad= zeros(21, ny, nx, num_of_layer+1);
lay_C_pow = zeros(21, ny, nx, num_of_layer+1);

%-- 1st: free surface
lay_elev   (:,:,1) =    -dep(1);
lay_den     (:,:,1) =     1e3;
lay_den_grad(:,:,1) =  den_grad(1);
lay_den_pow (:,:,1) =  den_pow (1);

lay_C      ( 1,:,:,1) = 25.2*1e9   ; %C11
lay_C      ( 3,:,:,1) = 10.9620*1e9; %C13
lay_C      (12,:,:,1) = 18.0*1e9   ; %C33
lay_C      (19,:,:,1) = 5.12*1e9   ; %C55
lay_C      (21,:,:,1) = 7.168*1e9  ; %C66
lay_C_grad (:,:,:,1) = C_grad(1);
lay_C_pow  (:,:,:,1) = C_pow (1);

%-- 2nd: basin
lay_den    (   :,:,2) =     1.5e3;
lay_den_grad(:,:,2) =  den_grad(2);
lay_den_pow (:,:,2) =  den_pow (2);
lay_C      ( 1,:,:,2) = 25.2*1e9    * 2; %C11
lay_C      ( 3,:,:,2) = 10.9620*1e9 * 2; %C13
lay_C      (12,:,:,2) = 18.0*1e9    * 2; %C33
lay_C      (19,:,:,2) = 5.12*1e9    * 2; %C55
lay_C      (21,:,:,2) = 7.168*1e9   * 2; %C66
lay_C_grad (:,:,:,2) = C_grad(2);
lay_C_pow  (:,:,:,2) = C_pow (2);

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
lay_den    (   :,:,3) =     1.2e3;
lay_den_grad(:,:,3) =  den_grad(3);
lay_den_pow (:,:,3) =  den_pow (3);
lay_C      ( 1,:,:,3) = 25.2*1e9    * 1.5; %C11
lay_C      ( 3,:,:,3) = 10.9620*1e9 * 1.5; %C13
lay_C      (12,:,:,3) = 18.0*1e9    * 1.5; %C33
lay_C      (19,:,:,3) = 5.12*1e9    * 1.5; %C55
lay_C      (21,:,:,3) = 7.168*1e9   * 1.5; %C66
lay_C_grad (:,:,:,3) = C_grad(3);
lay_C_pow  (:,:,:,3) = C_pow (3);

num_circle_x = 3; %- how many circle along x
num_circle_y = 1; %- how many circle along y

lay_elev(:,:,3) =    cos(num_circle_x * x2d / Lx * 2*pi)  ...
                  .* cos(num_circle_y * y2d/Ly*2*pi) ...
                  .* 5 * dx ...
                  - dep(3);

%-- 4rd: topo
lay_den    (   :,:,4) =     1.0e3;
lay_den_grad(:,:,4) =  den_grad(4);
lay_den_pow (:,:,4) =  den_pow (4);
lay_C      ( 1,:,:,4) = 25.2*1e9    ; %C11
lay_C      ( 3,:,:,4) = 10.9620*1e9 ; %C13
lay_C      (12,:,:,4) = 18.0*1e9    ; %C33
lay_C      (19,:,:,4) = 5.12*1e9    ; %C55
lay_C      (21,:,:,4) = 7.168*1e9   ; %C66
lay_C_grad (:,:,:,4) = C_grad(4);
lay_C_pow  (:,:,:,4) = C_pow (4);

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
media_type = 'elastic_aniso_cij'

fnm_ou = 'basin_el_aniso.md3lay'

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
              for n = 1 : 21
            	  fprintf(fid, ' %g', lay_C(n,j,i,ilay), C_grad(ilay), C_pow(ilay));
              end
              % return
            	fprintf(fid, '\n');
          end
      end
  end
fclose(fid);

%------------------------------------------------------------------------------
%-- create md3lay structure
%------------------------------------------------------------------------------

md.media_type = 'elastic_aniso_cij'
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

  md.C11     {n} = permute(squeeze(lay_C     ( 1,:,:,n)),[2,1]);
  md.C11_coef{n} = permute(squeeze(lay_C_grad( 1,:,:,n)),[2,1]);
  md.C11_pow {n} = permute(squeeze(lay_C_pow ( 1,:,:,n)),[2,1]);
  md.C12     {n} = permute(squeeze(lay_C     ( 2,:,:,n)),[2,1]);
  md.C12_coef{n} = permute(squeeze(lay_C_grad( 2,:,:,n)),[2,1]);
  md.C12_pow {n} = permute(squeeze(lay_C_pow ( 2,:,:,n)),[2,1]);
  md.C13     {n} = permute(squeeze(lay_C     ( 3,:,:,n)),[2,1]);
  md.C13_coef{n} = permute(squeeze(lay_C_grad( 3,:,:,n)),[2,1]);
  md.C13_pow {n} = permute(squeeze(lay_C_pow ( 3,:,:,n)),[2,1]);
  md.C14     {n} = permute(squeeze(lay_C     ( 4,:,:,n)),[2,1]);
  md.C14_coef{n} = permute(squeeze(lay_C_grad( 4,:,:,n)),[2,1]);
  md.C14_pow {n} = permute(squeeze(lay_C_pow ( 4,:,:,n)),[2,1]);
  md.C15     {n} = permute(squeeze(lay_C     ( 5,:,:,n)),[2,1]);
  md.C15_coef{n} = permute(squeeze(lay_C_grad( 5,:,:,n)),[2,1]);
  md.C15_pow {n} = permute(squeeze(lay_C_pow ( 5,:,:,n)),[2,1]);
  md.C16     {n} = permute(squeeze(lay_C     ( 6,:,:,n)),[2,1]);
  md.C16_coef{n} = permute(squeeze(lay_C_grad( 6,:,:,n)),[2,1]);
  md.C16_pow {n} = permute(squeeze(lay_C_pow ( 6,:,:,n)),[2,1]);
  md.C22     {n} = permute(squeeze(lay_C     ( 7,:,:,n)),[2,1]);
  md.C22_coef{n} = permute(squeeze(lay_C_grad( 7,:,:,n)),[2,1]);
  md.C22_pow {n} = permute(squeeze(lay_C_pow ( 7,:,:,n)),[2,1]);
  md.C23     {n} = permute(squeeze(lay_C     ( 8,:,:,n)),[2,1]);
  md.C23_coef{n} = permute(squeeze(lay_C_grad( 8,:,:,n)),[2,1]);
  md.C23_pow {n} = permute(squeeze(lay_C_pow ( 8,:,:,n)),[2,1]);
  md.C24     {n} = permute(squeeze(lay_C     ( 9,:,:,n)),[2,1]);
  md.C24_coef{n} = permute(squeeze(lay_C_grad( 9,:,:,n)),[2,1]);
  md.C24_pow {n} = permute(squeeze(lay_C_pow ( 9,:,:,n)),[2,1]);
  md.C25     {n} = permute(squeeze(lay_C     (10,:,:,n)),[2,1]);
  md.C25_coef{n} = permute(squeeze(lay_C_grad(10,:,:,n)),[2,1]);
  md.C25_pow {n} = permute(squeeze(lay_C_pow (10,:,:,n)),[2,1]);
  md.C26     {n} = permute(squeeze(lay_C     (11,:,:,n)),[2,1]);
  md.C26_coef{n} = permute(squeeze(lay_C_grad(11,:,:,n)),[2,1]);
  md.C26_pow {n} = permute(squeeze(lay_C_pow (11,:,:,n)),[2,1]);
  md.C33     {n} = permute(squeeze(lay_C     (12,:,:,n)),[2,1]);
  md.C33_coef{n} = permute(squeeze(lay_C_grad(12,:,:,n)),[2,1]);
  md.C33_pow {n} = permute(squeeze(lay_C_pow (12,:,:,n)),[2,1]);
  md.C34     {n} = permute(squeeze(lay_C     (13,:,:,n)),[2,1]);
  md.C34_coef{n} = permute(squeeze(lay_C_grad(13,:,:,n)),[2,1]);
  md.C34_pow {n} = permute(squeeze(lay_C_pow (13,:,:,n)),[2,1]);
  md.C35     {n} = permute(squeeze(lay_C     (14,:,:,n)),[2,1]);
  md.C35_coef{n} = permute(squeeze(lay_C_grad(14,:,:,n)),[2,1]);
  md.C35_pow {n} = permute(squeeze(lay_C_pow (14,:,:,n)),[2,1]);
  md.C36     {n} = permute(squeeze(lay_C     (15,:,:,n)),[2,1]);
  md.C36_coef{n} = permute(squeeze(lay_C_grad(15,:,:,n)),[2,1]);
  md.C36_pow {n} = permute(squeeze(lay_C_pow (15,:,:,n)),[2,1]);
  md.C44     {n} = permute(squeeze(lay_C     (16,:,:,n)),[2,1]);
  md.C44_coef{n} = permute(squeeze(lay_C_grad(16,:,:,n)),[2,1]);
  md.C44_pow {n} = permute(squeeze(lay_C_pow (16,:,:,n)),[2,1]);
  md.C45     {n} = permute(squeeze(lay_C     (17,:,:,n)),[2,1]);
  md.C45_coef{n} = permute(squeeze(lay_C_grad(17,:,:,n)),[2,1]);
  md.C45_pow {n} = permute(squeeze(lay_C_pow (17,:,:,n)),[2,1]);
  md.C46     {n} = permute(squeeze(lay_C     (18,:,:,n)),[2,1]);
  md.C46_coef{n} = permute(squeeze(lay_C_grad(18,:,:,n)),[2,1]);
  md.C46_pow {n} = permute(squeeze(lay_C_pow (18,:,:,n)),[2,1]);
  md.C55     {n} = permute(squeeze(lay_C     (19,:,:,n)),[2,1]);
  md.C55_coef{n} = permute(squeeze(lay_C_grad(19,:,:,n)),[2,1]);
  md.C55_pow {n} = permute(squeeze(lay_C_pow (19,:,:,n)),[2,1]);
  md.C56     {n} = permute(squeeze(lay_C     (20,:,:,n)),[2,1]);
  md.C56_coef{n} = permute(squeeze(lay_C_grad(20,:,:,n)),[2,1]);
  md.C56_pow {n} = permute(squeeze(lay_C_pow (20,:,:,n)),[2,1]);
  md.C66     {n} = permute(squeeze(lay_C     (21,:,:,n)),[2,1]);
  md.C66_coef{n} = permute(squeeze(lay_C_grad(21,:,:,n)),[2,1]);
  md.C66_pow {n} = permute(squeeze(lay_C_pow (21,:,:,n)),[2,1]);
end

%------------------------------------------------------------------------------
%-- export
%------------------------------------------------------------------------------

filename = 'basin_el_aniso.md3lay'

md3lay_export(filename, md);


