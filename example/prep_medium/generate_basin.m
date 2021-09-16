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

%- roughly depth of each interface
dep     = [-1e3, basin_r0, 3e3, 10e3];
%dep     = [0, basin_r0, 3e3, 10e3];

Vp      = [ 1500, 3000, 5000, 8000];
Vp_grad = [ 0.2 ,  0.0, 0.0,  0 ];
Vp_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

Vs       = [ 1000, 1500, 2500, 3000 ];
Vs_grad  = [ 0.2 ,  0.0, 0.0,  0 ];
Vs_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

den      = [ 1000, 2000, 3000, 3500 ];
den_grad = [ 0.2 ,  0.0, 0.0,  0 ];
den_pow  = [ 1.0 ,  1.0, 1.0,  1 ];

[x2d, y2d] = meshgrid(x1d,y1d);

%-- calculate elevation of each interface
elev_int = zeros(ny, nx, num_of_layer+1);

%-- 1st: free surface
elev_int(:,:,1) = -dep(1);

%-- 2nd: basin
elev_int(:,:,2) = elev_int(:,:,1);
for j = 1 : ny
for i = 1 : nx
    r2d = sqrt((x1d(i) - basin_x0)^2 + (y1d(j) - basin_y0)^2);
    %if in basin
    if r2d <= basin_r0
      %-- (x-x0)^2 + (y-y0)^2 + (e-e0)^2 = r^2
      elev_int(j,i,2) = - (sqrt( basin_r0^2 - r2d^2 ) + basin_e0 );
    end
end
end

%-- 3rd: topo
num_circle_x = 3; %- how many circle along x
num_circle_y = 1; %- how many circle along y

elev_int(:,:,3) =    cos(num_circle_x * x2d / Lx * 2*pi)  ...
                  .* cos(num_circle_y * y2d/Ly*2*pi) ...
                  .* 5 * dx ...
                  - dep(3);

%-- 4rd: topo
elev_int(:,:,4) = -dep(4);

%------------------------------------------------------------------------------
%-- plot
%------------------------------------------------------------------------------

figure;
for ilay = 1 : num_of_layer+1
    mesh(x2d, y2d, elev_int(:,:,ilay));
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

fnm_ou = 'basin_rho.md3lay'
fid = fopen(fnm_ou,'w');
fprintf(fid, '%d\n', num_of_layer + 1);
fprintf(fid, '%d %d %f %f %f %f\n', nx, ny, x0, y0, dx, dy);
  for ilay = 1 : num_of_layer+1
      for j = 1 : ny
          for i = 1 : nx
            	fprintf(fid, '%f %f %f %f\n', elev_int(j,i,ilay), ...
                    den(ilay), den_grad(ilay), den_pow(ilay));
          end
      end
  end
fclose(fid);

fnm_ou = 'basin_Vp.md3lay'
fid = fopen(fnm_ou,'w');
fprintf(fid, '%d\n', num_of_layer + 1);
fprintf(fid, '%d %d %f %f %f %f\n', nx, ny, x0, y0, dx, dy);
  for ilay = 1 : num_of_layer+1
      for j = 1 : ny
          for i = 1 : nx
            	fprintf(fid, '%f %f %f %f\n', elev_int(j,i,ilay), ...
                    Vp(ilay), Vp_grad(ilay), Vp_pow(ilay));
          end
      end
  end
fclose(fid);

fnm_ou = 'basin_Vs.md3lay'
fid = fopen(fnm_ou,'w');
fprintf(fid, '%d\n', num_of_layer + 1);
fprintf(fid, '%d %d %f %f %f %f\n', nx, ny, x0, y0, dx, dy);
  for ilay = 1 : num_of_layer+1
      for j = 1 : ny
          for i = 1 : nx
            	fprintf(fid, '%f %f %f %f\n', elev_int(j,i,ilay), ...
                    Vs(ilay), Vs_grad(ilay), Vs_pow(ilay));
          end
      end
  end
fclose(fid);

