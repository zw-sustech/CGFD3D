clear all;

draw_out_model = 1;

filename = 'can4_rho.md3lay'
fid = fopen(filename, 'r');

NI = fscanf(fid, '%d\n', 1); 
NX = fscanf(fid, '%d\n', 1); 
NY = fscanf(fid, '%d\n', 1); 
MINX = fscanf(fid, '%f\n', 1); 
MINY = fscanf(fid, '%f\n', 1); 
DX   = fscanf(fid, '%f\n', 1); 
DY   = fscanf(fid, '%f\n', 1); 
 
for ni = 1:NI
    for j = 1:NY
        for i = 1:NX
            elevation(j,i,ni) = fscanf(fid, '%f',1);
            data(j,i,ni) = fscanf(fid, '%f', 1);
            par_grad(j,i,ni) = fscanf(fid, '%f', 1);
            par_pow(j,i, ni) = fscanf(fid, '%f', 1);
        end
    end
end

fclose(fid);

xvec = [0:NX-1]*DX + MINX;
yvec = [0:NY-1]*DY + MINY;
[X,Y] = meshgrid(xvec,yvec);
% plot
figure;
for ni = 1:NI
    mesh(X, Y, elevation(:,:,ni));
    hold on;
    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    title('The given model');
end

if draw_out_model
    nx = 501;  ny = 501; nz = 450;
    x0 = -2000; y0 = 0; z0 = -4500;
    dx = 10; dy = 10; dz = 10;
    data_file = 'rho.dat';
 %   xslice = [-2000,0,2000];
 xslice = [];
    yslice = [5000];
    zslice = [];
    drawmodel(nx, ny, nz, x0, y0, z0, dx, dy, dz, data_file, xslice, yslice, zslice);
end
