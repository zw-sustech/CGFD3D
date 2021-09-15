%function drawmodel(nx,ny,nz,x0,y0,z0,dx, dy,dz, data_file, xslice, yslice, zslice)
close all;
clear all;
    nx = 501;  ny = 501; nz = 450;
    x0 = -2000; y0 = 0; z0 = -4500;
    dx = 10; dy = 10; dz = 10;

 %   xslice = [-2000,0,2000];
    xslice = [];
    yslice = [1000];
    zslice = [];

xvec = x0+[0:nx-1]*dx;
yvec = y0+[0:ny-1]*dy;
zvec = z0+[0:nz-1]*dz;

data_file1 = 'rho1.dat';
data_file3 = 'rho3.dat';
data_file2 = 'rho2.dat';

%=== kappa, mu har
fid = fopen(data_file1);
data1 = fread(fid,nx*ny*nz,'float');
fclose(fid);

%fid = fopen(data_file3);
%data3 = fread(fid,nx*ny*nz,'float');
%fclose(fid);

%fid = fopen(data_file2);
%data2 = fread(fid,nx*ny*nz,'float');
%fclose(fid);

u1 = permute(reshape(data1, nx, ny, nz),[2 1 3]);
%u3 = permute(reshape(data3, nx, ny, nz),[2 1 3]);
%u2 = permute(reshape(data2, nx, ny, nz),[2 1 3]);

%figure;
figure;
h = slice(xvec,yvec,zvec,u1,xslice,yslice,zslice);
set(h, 'edgecolor','none');
grid off;
hold on;

%colorbar;
j = 1;
if 00
    for i = 1:nx*ny*nz
        if (data1(i) ~= data2(i))
            k(j) = i
            j = j+1;
        end
    end
end

%end