% Generate a modified Can4 model (from E2VP)
% reference: Chaljub et al, 2015, 3D numerical simulation of earthquake ground motion
%   in sedimentary basins: testing accuracy through stringent models
clear all;
close all;

if_write_rho = 1;
if_write_vp = 1;
if_write_vs = 1;

NI = 4;

H(1) =     0;
H(2) = - 173*3;
H(3) = - 725*3;
H(4) = -1156*3;

% media info od the first layer
rho(1) = 1000;
vp(1)  = 1500;
vs(1)  = 1000;
par_grad = 0.0;
par_pow  = 1; 


for ni = 2:NI
    rho(ni) = rho(ni-1) * 1.5;
    vp(ni)  =  vp(ni-1) * 1.5;
    vs(ni)  =  vs(ni-1) * 1.5;
end


%MINX = -2500.0; DX = 500.0; NX = 11;
MINX = -2500.0; DX = 500.0; NX = 101;
MINY = -2500.0; DY = 500.0; NY = 101;

xvec = [0:NX-1]*DX + MINX;
yvec = [0:NY-1]*DY + MINY;
[X,Y] = meshgrid(xvec,yvec);

% calculate elevation 
for ni = 1:NI
    for j =1:NY
        for i = 1:NX
            elevation(j,i,ni) = H(ni) * min(1, (10000-xvec(i))/5000);
            if (elevation(j,i,ni) > 0) 
               elevation(j,i,ni) = 0;
            end
        end
    end
end

% plot
figure;
for ni = 1:NI
    mesh(X, Y, elevation(:,:,ni));
    hold on;
    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    title('The build model');
end

%%%% 
%figure;
%drawmodel(501,501,300,-2500,0,-3000,10,10,10,'rho.dat',[],0,[]);

if if_write_rho
	fid = fopen('can4_rho.md3lay','w');
	fprintf(fid, '%d\n', NI);
	fprintf(fid, '%d %d %f %f %f %f\n', NX, NY, MINX, MINY, DX, DY);
    for ni = 1:NI
        for j = 1:NY
            for i = 1:NX
              	fprintf(fid, '%f %f %f %f\n', elevation(j,i,ni), rho(ni), par_grad, par_pow);
            end
        end
    end
    fclose(fid);
end

if if_write_vp
	fid = fopen('can4_vp.md3lay','w');
	fprintf(fid, '%d\n', NI);
	fprintf(fid, '%d %d %f %f %f %f\n', NX, NY, MINX, MINY, DX, DY);
    for ni = 1:NI
        for j = 1:NY
            for i = 1:NX
              	fprintf(fid, '%f %f %f %f\n', elevation(j,i,ni), vp(ni), par_grad, par_pow);
            end
        end
    end
    fclose(fid);
end

if if_write_vs
	fid = fopen('can4_vs.md3lay','w');
	fprintf(fid, '%d\n', NI);
	fprintf(fid, '%d %d %f %f %f %f\n', NX, NY, MINX, MINY, DX, DY);
    for ni = 1:NI
        for j = 1:NY
            for i = 1:NX
              	fprintf(fid, '%f %f %f %f\n', elevation(j,i,ni), vs(ni), par_grad, par_pow);
            end
        end
    end
    fclose(fid);
end
