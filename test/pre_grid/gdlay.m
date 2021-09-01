%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Example of generating the interface file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc

%%
layer_nx = 100;
layer_ny = 100;
layer_number = 4; % Number of interfaces 
NCellPerlay = [ 25 10 25 ]; % Number of nodes per layer 
VmapSpacingIsequal =  [1 0 1]; % 1:The grid spacing is equal; 0: Is not.

num_grid = layer_nx*layer_ny*layer_number;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% generate interface grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx = 1:layer_nx;
yy = 1:layer_ny;
[GridY,GridX] = meshgrid(yy,xx);
layer_depth = [-6000 -2500 -1500 0];

%%%% hill
[m1,m2]=meshgrid(yy-layer_ny/2, xx-layer_nx/2);
r=sqrt( m1.^2+m2.^2 );
z0 =exp(-r.^2/(layer_nx/4)^2);

mx = zeros(layer_nx,layer_ny,layer_number);
my = zeros(layer_nx,layer_ny,layer_number);
mz = zeros(layer_nx,layer_ny,layer_number);

for i = 1:layer_number
    mx(:,:,i) = GridX*100;
    my(:,:,i) = GridY*100;
    mz(:,:,i) = z0*5000*((i-1)/layer_number) + layer_depth(i);
end

layer3d = zeros(num_grid*3,1);
nn = 1;
for i=1:layer_number
    for j=1:layer_ny
        for k = 1:layer_nx
            layer3d(nn           ) = mx(k,j,i);
            layer3d(nn+num_grid  ) = my(k,j,i);
            layer3d(nn+num_grid*2) = mz(k,j,i);
            nn = nn+1;
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Write gdlay text  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen("test_grid.gdlay","w"); % Output file name 
fprintf(fid,'%6d\n',layer_number);
fprintf(fid,'%6d %6d %6d\n',NCellPerlay);
fprintf(fid,'%6d %6d %6d\n',VmapSpacingIsequal);
fprintf(fid,'%6d %6d\n',layer_nx, layer_ny);

for i=1: num_grid
    fprintf(fid,'%12.2f %12.2f %12.2f\n',[layer3d(i) layer3d(i+num_grid) layer3d(i+num_grid*2)]);
end
fclose(fid);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
hold off
for ii = 1:size(mx,3)
    surf(mx(:,:,ii),my(:,:,ii),mz(:,:,ii))
    hold on
end
colorbar();
shading interp
% set(gca,'zdir','reverse')
set(gca,'fontsize',16)
ylabel('Y');
xlabel('X');
zlabel('Depth','fontsize',20);
set(gcf,'Position',[50 50 900 1000])
axis equal tight
box on


