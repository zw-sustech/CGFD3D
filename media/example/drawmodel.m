function drawmodel(nx,ny,nz,x0,y0,z0,dx, dy,dz, data_file, xslice, yslice, zslice)

xvec = x0+[0:nx-1]*dx;
yvec = y0+[0:ny-1]*dy;
zvec = z0+[0:nz-1]*dz;

%=== kappa, mu har
fid = fopen(data_file);
data = fread(fid,nx*ny*nz,'float');
fclose(fid);

u = permute(reshape(data, nx, ny, nz),[2 1 3]);

%figure;
h = slice(xvec,yvec,zvec,u,xslice,yslice,zslice);
set(h, 'edgecolor','none');
grid off;
hold on;

colorbar;

end