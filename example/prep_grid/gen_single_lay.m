%------------------------------------------------------------------------------
%-- Example of generating the interface file 
%------------------------------------------------------------------------------

clear

%------------------------------------------------------------------------------
%-- generate x and y
%------------------------------------------------------------------------------

%-- should 6 more than FD points due to ghosts points
nghost = 3;
nx = 326;
ny = 306;
x0 = 0.0;
y0 = 0.0;
dx = 100.0;
dy = 100.0;

x1d = [0 : nx-1] * dx + x0 - nghost * dx;
y1d = [0 : ny-1] * dy + y0 - nghost * dy;

flag_figure = 1;

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

if flag_figure == 1

  figure
  surf(y1d,x1d,Z,'edgecolor','none');
  xlabel('y-axis');
  ylabel('x-axis');
  zlabel('z-axis');
  daspect([1,1,1]);
  camlight;
  title('generated topography');
  print(gcf,'-dpng','gdlay_topography.png');
end

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

num_of_interfaces = 2; 

% 1 less than num_of_interfaces, total cell should be nghost more than FD points
num_of_cell_per_layer = [ 63 ];
dz_is_equal_of_layer  = [ 1 ]; % 1:The grid spacing is equal; 0: Is not.
avg_dz_of_layer       = [ 100 ];
smo_sigma             = [ 10, 0]; % last one is not used

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
  topo  = imgaussfilt(free_topo, smo_sigma(ilay), ...
                      'FilterSize',2*smo_width(ilay)+1);
  z3d(:,:,ilay) = topo + z_of_interfaces(ilay);
end

%------------------------------------------------------------------------------
%-- plot
%------------------------------------------------------------------------------

if flag_figure == 1

  figure
  for ilay = 1 : num_of_interfaces
    surf(x3d(:,:,ilay), y3d(:,:,ilay), z3d(:,:,ilay))
    hold on
  end
  xlabel('y-axis');
  ylabel('x-axis');
  zlabel('z-axis');
  shading flat
  colorbar
  camlight
  %daspect([8,10,110e3*8/5]);
  daspect([1,1,1]);
  title('layers in gdlay');
  print(gcf,'-dpng','gdlay_layers.png');

  for ilay = 1 : num_of_interfaces
    figure
    surf(x3d(:,:,ilay), y3d(:,:,ilay), z3d(:,:,ilay))
    xlabel('y-axis');
    ylabel('x-axis');
    zlabel('z-axis');
    shading flat
    colorbar
    camlight
    %daspect([8,10,110e3*8/5]);
    daspect([1,1,1]);
    title(['layer ',num2str(ilay),'th']);
    print(gcf,'-dpng',['gdlay_lay',num2str(ilay),'.png']);
  end

end

%==============================================================================
%-- write .gdlay file
%==============================================================================
gdlay_file = 'random_topo_single.gdlay'

gdlay_export(gdlay_file, ...
             num_of_interfaces, nx, ny, ...
             num_of_cell_per_layer, ...
             dz_is_equal_of_layer, ...
             x3d, y3d, z3d);

