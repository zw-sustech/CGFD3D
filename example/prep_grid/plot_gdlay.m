%------------------------------------------------------------------------------
%-- Example of generating the interface file 
%------------------------------------------------------------------------------

clear

%------------------------------------------------------------------------------
%-- read .gdlay file
%------------------------------------------------------------------------------

gdlay_file = 'random_topo.gdlay'
%gdlay_file = 'random_topo_single.gdlay'

gd = gdlay_import(gdlay_file);

%------------------------------------------------------------------------------
%-- plot surface
%------------------------------------------------------------------------------

figure
hold off
for n = 1 : gd.num_of_interfaces
    surf(gd.x3d(:,:,n), gd.y3d(:,:,n), gd.z3d(:,:,n));
    %plot3(gd.x3d(:,:,n), gd.y3d(:,:,n), gd.z3d(:,:,n));
    hold on
end

colorbar();
shading flat
% set(gca,'zdir','reverse')
set(gca,'fontsize',16)
title(['grid layers of ', gdlay_file]);
ylabel('X');
xlabel('Y');
zlabel('Depth','fontsize',20);
set(gcf,'Position',[50 50 900 1000])
axis equal tight
box on

%------------------------------------------------------------------------------
%-- plot slice
%------------------------------------------------------------------------------

i_slice = nearest( gd.nx / 2 );
j_slice = nearest( gd.ny / 2 );

figure
plot3(squeeze(gd.x3d(i_slice,:,:)), ...
     squeeze(gd.y3d(i_slice,:,:)), ...
     squeeze(gd.z3d(i_slice,:,:)), 'k');

hold on
plot3(squeeze(gd.x3d(:,j_slice,:)), ...
     squeeze(gd.y3d(:,j_slice,:)), ...
     squeeze(gd.z3d(:,j_slice,:)), 'b');


colorbar();
shading flat
% set(gca,'zdir','reverse')
set(gca,'fontsize',16)
title(['grid layers throught two vertical slices']);
ylabel('X');
xlabel('Y');
zlabel('Depth','fontsize',20);
set(gcf,'Position',[50 50 900 1000])
axis equal tight
box on

%------------------------------------------------------------------------------
%-- create demo grid and plot grid on two vertical slices
%------------------------------------------------------------------------------

grid_nx = gd.nx;
grid_ny = gd.ny;
grid_nz = sum(gd.num_of_cell_per_layer) + 1;

grid_x3d = zeros(grid_nx, grid_ny, grid_nz);
grid_y3d = zeros(grid_nx, grid_ny, grid_nz);
grid_z3d = zeros(grid_nx, grid_ny, grid_nz);

%-- loop to set points of each layer
for n = 1 : gd.num_of_interfaces - 1
  % start k of this layer
  if n == 1
   k1 = 1;
  else
   k1 = sum(gd.num_of_cell_per_layer(1:n-1)) + 1;
  end
  % end k of this layer
  k2 = sum(gd.num_of_cell_per_layer(1:n)) + 1;
  
  %-- start layer
  grid_x3d(:,:,k1) = gd.x3d(:,:,n);
  grid_y3d(:,:,k1) = gd.y3d(:,:,n);
  grid_z3d(:,:,k1) = gd.z3d(:,:,n);

  %-- end layer
  grid_x3d(:,:,k2) = gd.x3d(:,:,n+1);
  grid_y3d(:,:,k2) = gd.y3d(:,:,n+1);
  grid_z3d(:,:,k2) = gd.z3d(:,:,n+1);

  %-- interp inner point
  for k_cell = 1 : gd.num_of_cell_per_layer(n)
    % interp coef
    L2 = (k_cell - 1) / gd.num_of_cell_per_layer(n);
    L1 = 1.0 - L2;
    % index
    k = k1 + k_cell - 1;

    grid_x3d(:,:,k) = grid_x3d(:,:,k1) * L1 + grid_x3d(:,:,k2) * L2;
    grid_y3d(:,:,k) = grid_y3d(:,:,k1) * L1 + grid_y3d(:,:,k2) * L2;
    grid_z3d(:,:,k) = grid_z3d(:,:,k1) * L1 + grid_z3d(:,:,k2) * L2;
  end
end

i_slice = nearest( gd.nx / 2 );
j_slice = nearest( gd.ny / 2 );

%-- plot slice grid
figure
surf(squeeze(grid_x3d(i_slice,:,:)), ...
     squeeze(grid_y3d(i_slice,:,:)), ...
     squeeze(grid_z3d(i_slice,:,:)));

hold on
surf(squeeze(grid_x3d(:,j_slice,:)), ...
     squeeze(grid_y3d(:,j_slice,:)), ...
     squeeze(grid_z3d(:,j_slice,:)));
shading faceted;

colorbar();
shading faceted
% set(gca,'zdir','reverse')
set(gca,'fontsize',16)
title('equal interpolated grid along two vertical slices');
ylabel('X');
xlabel('Y');
zlabel('Depth','fontsize',20);
set(gcf,'Position',[50 50 900 1000])
axis equal tight
box on

%-- plot slice grid step
grid_dz = grid_z3d(:,:,2:grid_nz) - grid_z3d(:,:,1:grid_nz-1);

figure
surf(squeeze(grid_x3d(i_slice,:,1:end-1)), ...
     squeeze(grid_y3d(i_slice,:,1:end-1)), ...
     squeeze(grid_z3d(i_slice,:,1:end-1)), ...
     squeeze(grid_dz (i_slice,:,:)));

hold on
surf(squeeze(grid_x3d(:,j_slice,1:end-1)), ...
     squeeze(grid_y3d(:,j_slice,1:end-1)), ...
     squeeze(grid_z3d(:,j_slice,1:end-1)), ...
     squeeze(grid_dz (:,j_slice,:)));
shading faceted;

colorbar();
shading faceted
% set(gca,'zdir','reverse')
set(gca,'fontsize',16)
title('averaged vertical grid spacing along two vertical slices');
ylabel('X');
xlabel('Y');
zlabel('Depth','fontsize',20);
set(gcf,'Position',[50 50 900 1000])
axis equal tight
box on
