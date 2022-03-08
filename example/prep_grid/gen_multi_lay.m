%------------------------------------------------------------------------------
%-- Example of generating the interface file 
%------------------------------------------------------------------------------

clear

%------------------------------------------------------------------------------
%-- generate x and y
%------------------------------------------------------------------------------

%-- should 6 more than FD points due to ghosts points
nghost = 3;
nx = 126;
ny = 106;
x0 = 0.0;
y0 = 0.0;
dx = 100.0;
dy = 100.0;

x1d = [0 : nx-1] * dx + x0 - nghost * dx;
y1d = [0 : ny-1] * dy + y0 - nghost * dy;

%------------------------------------------------------------------------------
%-- generate or load topography
%------------------------------------------------------------------------------

topo = zeros(nx,ny);

%-- from https://www.mathworks.com/matlabcentral/answers/218806-random-gaussian-surface-generation
N = [nx ny]; % size in pixels of image
F = 10;        % frequency-filter width
[X,Y] = ndgrid(1:N(1),1:N(2));
i = min(X-1,N(1)-X+1);
j = min(Y-1,N(2)-Y+1);
H = exp(-.5*(i.^2+j.^2)/F^2);
Z = real(ifft2(H.*fft2(randn(N))));
Z = Z .* 1e3;
figure
surf(X,Y,Z,'edgecolor','none');
light;

%-- get from above surface
%for j = 1 : ny
%for i = 1 : nx
%  free_topo(i,j) = Z(i,j);
%end
%end
free_topo = Z;

%------------------------------------------------------------------------------
%-- generate grid interfaces
%------------------------------------------------------------------------------

%-- note: from bottom to top

num_of_interfaces = 4; 

% 1 less than num_of_interfaces, total cell should be nghost more than FD points
num_of_cell_per_layer = [ 43 10 10 ];
dz_is_equal_of_layer  = [ 1 0 1]; % 1:The grid spacing is equal; 0: Is not.
avg_dz_of_layer       = [ 100, 75, 50 ];
smooth_length         = [ 20, 5, 1, 1];

%-- use avg_dz and num_of_cell to esti z-axis of each interface
%-- last elem is the free surface
z_of_interfaces(num_of_interfaces) = 0;
for ilay = num_of_interfaces-1 : -1 : 1
  z_of_interfaces(ilay) = z_of_interfaces(ilay+1) ...
        - avg_dz_of_layer(ilay) * num_of_cell_per_layer(ilay);
end

%-- construct grid interfaces from free_topo and z_of_interfaces
x3d = zeros(nx,ny,num_of_interfaces);
y3d = zeros(nx,ny,num_of_interfaces);
z3d = zeros(nx,ny,num_of_interfaces);

%-- set x3d
for n = 1 : num_of_interfaces
for j = 1 : ny
  x3d(:,j,n) = x1d;
end
end

%-- set y3d
for n = 1 : num_of_interfaces
for i = 1 : nx
  y3d(i,:,n) = y1d';
end
end

%- first same to free surface
z3d(:,:,num_of_interfaces) = free_topo;

for ilay = num_of_interfaces-1 : -1 : 1
  %-- smooth free_topo 
  topo  = smooth2(free_topo, smooth_length(ilay));
  z3d(:,:,ilay) = topo + z_of_interfaces(ilay);
end

%------------------------------------------------------------------------------
%-- plot
%------------------------------------------------------------------------------

%figure
%hold off
%for n = 1 : num_of_interfaces
%    %surf(x1d, y1d, permute(z3d(:,:,n),[2 1]));
%    surf(x3d(:,:,n), y3d(:,:,n), z3d(:,:,n));
%    hold on
%end
%
%colorbar();
%shading interp
%% set(gca,'zdir','reverse')
%set(gca,'fontsize',16)
%ylabel('X');
%xlabel('Y');
%zlabel('Depth','fontsize',20);
%set(gcf,'Position',[50 50 900 1000])
%axis equal tight
%box on

%==============================================================================
%-- write .gdlay file
%==============================================================================
gdlay_file = 'random_topo.gdlay'

gdlay_export(gdlay_file, ...
             num_of_interfaces, nx, ny, ...
             num_of_cell_per_layer, ...
             dz_is_equal_of_layer, ...
             x3d, y3d, z3d);

