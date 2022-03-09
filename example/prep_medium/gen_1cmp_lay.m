% Generate a modified Can4 model (from E2VP)
% reference: Chaljub et al, 2015, 3D numerical simulation of earthquake ground motion
%   in sedimentary basins: testing accuracy through stringent models
clear all;
close all;

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

elevation = zeros(NY,NX,NI);
grad3d    = zeros(NY,NX,NI);
pow3d     = zeros(NY,NX,NI);

Vp3d     = zeros(NY,NX,NI);
Vs3d     = zeros(NY,NX,NI);
Dp3d     = zeros(NY,NX,NI);

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

    % layer values
    Vp3d(:,:,ni) = vp(ni);
    Vs3d(:,:,ni) = vs(ni);
    Dp3d(:,:,ni) = rho(ni);
end
% const grad and pow
grad3d(:,:,:) = par_grad;
pow3d(:,:,:) = par_pow;

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

%------------------------------------------------------------------------------
%-- create md3lay structure
%------------------------------------------------------------------------------

md.media_type = 'one_component'
%  one_component, 
%  acoustic_isotropic, 
%  elastic_isotropic, 
%  elastic_vti_prem, elastic_vti_thomsen, elastic_vti_cij,
%  elastic_tti_thomsen, elastic_tti_bond,
%  elastic_aniso_cij

md.num_of_intfce = NI;

md.nx = NX;
md.ny = NY;
md.dx = DX;
md.dy = DY;
md.x0 = MINX;
md.y0 = MINY;

%-- permute if order is not [x,y,layer]
for n = 1 : md.num_of_intfce
  md.elev{n}    = permute(elevation(:,:,n),[2,1]);

  md.val     {n} = permute(Dp3d     (:,:,n),[2,1]);
  md.val_coef{n} = permute(grad3d(:,:,n),[2,1]);
  md.val_pow {n} = permute(pow3d (:,:,n),[2,1]);
end

%-- export rho
filename = 'can4_rho.md3lay'
md3lay_export(filename, md);

%-- export Vp
for n = 1 : md.num_of_intfce
  md.val{n} = permute(Vp3d(:,:,n),[2,1]);
end
filename = 'can4_vp.md3lay'
md3lay_export(filename, md);

%-- export Vs
for n = 1 : md.num_of_intfce
  md.val{n} = permute(Vs3d(:,:,n),[2,1]);
end
filename = 'can4_vs.md3lay'
md3lay_export(filename, md);

