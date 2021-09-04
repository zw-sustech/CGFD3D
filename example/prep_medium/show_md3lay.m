clear all;

filename = 'basin_Vp.md3lay'

%-- read in data
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

%-- plot

figure;
for ni = 1:NI
    mesh(X, Y, elevation(:,:,ni));
    hold on;
    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    title('The given model');
end

