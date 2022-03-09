clear all;

fnm = 'topolay_el_iso.md3grd'

%-- read file
md = md3grd_import(fnm);

xvec = [0:md.nx-1]* md.dx + md.x0;
yvec = [0:md.ny-1]* md.dy + md.y0;

[x2d, y2d] = meshgrid(xvec,yvec);

% creat structure mesh
for n = 1 : md.num_of_intfce
  md.x3d{n} = zeros(md.ny, md.nx, md.num_of_point_per_lay(n));
  md.y3d{n} = zeros(md.ny, md.nx, md.num_of_point_per_lay(n));
  for k = 1 : md.num_of_point_per_lay(n)
    md.x3d{n}(:,:,k) = x2d;
    md.y3d{n}(:,:,k) = y2d;
  end
end

%------------------------------------------------------------------------------
%-- plot slice
%------------------------------------------------------------------------------

i_slice = nearest( md.nx / 2 );
j_slice = nearest( md.ny / 2 );

figure;

for n = 1 : md.num_of_intfce
    surf( squeeze( md.x3d{n}(j_slice,:,:)), ...
          squeeze( md.y3d{n}(j_slice,:,:)), ...
          squeeze(md.elev{n}(:,j_slice,:)), ...
          squeeze(md.Vp  {n}(:,j_slice,:)));
    hold on;

    surf( squeeze( md.x3d{n}(:,i_slice,:)), ...
          squeeze( md.y3d{n}(:,i_slice,:)), ...
          squeeze(md.elev{n}(i_slice,:,:)), ...
          squeeze(md.Vp  {n}(i_slice,:,:)));

    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    shading faceted;
    %shading flat;
    %title('structure grid');
    colorbar;
    title('Vp');
end

%------------------------------------------------------------------------------
%-- plot layers
%------------------------------------------------------------------------------

figure;

for n = 1 : md.num_of_intfce
    surf( squeeze( md.x3d{n}(:,:,1)), ...
          squeeze( md.y3d{n}(:,:,1)), ...
          squeeze(permute(md.elev{n}(:,:,1),[2 1])), ...
          squeeze(permute(md.Vp  {n}(:,:,1),[2 1])) );
    hold on;

    xlabel('x','fontsize', 12);
    ylabel('y','fontsize', 12);
    axis image;
    %shading flat;
    shading faceted;
    %title('structure grid');
    colorbar;
    title('Vp at top of each layers');
end
