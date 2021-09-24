clear all;

filename = 'basin_el_iso.md3lay'

%-- read in data
fid = fopen(filename, 'r');

media_type = fscanf(fid, '%s\n', 1); 
NI = fscanf(fid, '%d\n', 1); 
NX = fscanf(fid, '%d\n', 1); 
NY = fscanf(fid, '%d\n', 1); 
MINX = fscanf(fid, '%f\n', 1); 
MINY = fscanf(fid, '%f\n', 1); 
DX   = fscanf(fid, '%f\n', 1); 
DY   = fscanf(fid, '%f\n', 1); 

lay_elev    = zeros(NY, NX, NI);
lay_den     = zeros(NY, NX, NI);
lay_den_grad= zeros(NY, NX, NI);
lay_den_pow = zeros(NY, NX, NI);
lay_Vp      = zeros(NY, NX, NI);
lay_Vp_grad = zeros(NY, NX, NI);
lay_Vp_pow  = zeros(NY, NX, NI);
lay_Vs      = zeros(NY, NX, NI);
lay_Vs_grad = zeros(NY, NX, NI);
lay_Vs_pow  = zeros(NY, NX, NI);


if strcmp(media_type, 'acoustic_isotropic') == 1
   for ni = 1:NI
       for j = 1:NY
           for i = 1:NX
               lay_elev    (j,i,ni) = fscanf(fid, '%f',1);
               lay_den     (j,i,ni) = fscanf(fid, '%f', 1);
               lay_den_grad(j,i,ni) = fscanf(fid, '%f', 1);
               lay_den_pow (j,i,ni) = fscanf(fid, '%f', 1);
               lay_Vp     (j,i,ni) = fscanf(fid, '%f', 1);
               lay_Vp_grad(j,i,ni) = fscanf(fid, '%f', 1);
               lay_Vp_pow (j,i,ni) = fscanf(fid, '%f', 1);
           end
       end
   end
elseif strcmp(media_type, 'elastic_isotropic') == 1
   for ni = 1:NI
       for j = 1:NY
           for i = 1:NX
               lay_elev    (j,i,ni) = fscanf(fid, '%f',1);
               lay_den     (j,i,ni) = fscanf(fid, '%f', 1);
               lay_den_grad(j,i,ni) = fscanf(fid, '%f', 1);
               lay_den_pow (j,i,ni) = fscanf(fid, '%f', 1);
               lay_Vp     (j,i,ni) = fscanf(fid, '%f', 1);
               lay_Vp_grad(j,i,ni) = fscanf(fid, '%f', 1);
               lay_Vp_pow (j,i,ni) = fscanf(fid, '%f', 1);
               lay_Vs     (j,i,ni) = fscanf(fid, '%f', 1);
               lay_Vs_grad(j,i,ni) = fscanf(fid, '%f', 1);
               lay_Vs_pow (j,i,ni) = fscanf(fid, '%f', 1);
           end
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
    mesh(X, Y, lay_elev(:,:,ni));
    hold on;
    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    title('elevation');
end

figure;
for ni = 1:NI
    surf(X, Y, lay_elev(:,:,ni), lay_Vp(:,:,ni));
    hold on;
    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    shading flat;
    title('Vp');
end
