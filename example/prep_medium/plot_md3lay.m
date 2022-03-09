clear all;

%filename = 'basin_el_iso.md3lay'
%filename = 'basin_ac_iso.md3lay'
%filename = 'basin_el_aniso.md3lay'
%filename = 'basin_el_tti.md3lay'
%filename = 'basin_el_vti.md3lay'
filename = 'can4_vp.md3lay'

%-- read file
md = md3lay_import(filename);

xvec = [0:md.nx-1]* md.dx + md.x0;
yvec = [0:md.ny-1]* md.dy + md.y0;

[x2d, y2d] = meshgrid(xvec,yvec);

%------------------------------------------------------------------------------
%-- plot slice
%------------------------------------------------------------------------------

figure;

for n = 1 : md.num_of_intfce
    mesh(x2d, y2d, squeeze(permute(md.elev{n}(:,:),[2,1])));
    hold on;
    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    title('elevation of each interfaces');
end

figure;

for n = 1 : md.num_of_intfce
    surf(x2d, y2d, ...
         squeeze(permute(md.elev{n}(:,:),[2,1])),  ...
         squeeze(permute(md.val  {n}(:,:),[2,1]))  ...
         );
    hold on;
    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    title('Vp of each interfaces');
end

%for n = 1 : md.num_of_intfce
%    surf(x2d, y2d, ...
%         squeeze(permute(md.elev{n}(:,:),[2,1])),  ...
%         squeeze(permute(md.C11 {n}(:,:),[2,1]))  ...
%         );
%    hold on;
%    xlabel('x','fontsize', 12);
%    ylabel('y','fontsize', 12);
%    axis image;
%    title('C11 of each interfaces');
%end
