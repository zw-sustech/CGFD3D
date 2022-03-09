function md = md3grd_import(fnm_in)

% This function exports structure model to md3grd file.
%  fnm_in: string, the md3grd file name;
%  md: model structure. What should be kept in md depends on its media_type:
%     media_type: supported value 
%       one_component, 
%       acoustic_isotropic, 
%       elastic_isotropic, 
%       elastic_vti_prem, elastic_vti_thomsen, elastic_vti_cij,
%       elastic_tti_thomsen, elastic_tti_bond,
%       elastic_aniso_cij
%
%     num_of_intfce: number of interfaces or layers,
%     nx: number of sampling points along x-axis
%     ny: number of sampling points along y-axis
%     x0: x0 of first sampling points along x-axis
%     y0: y0 of first sampling points along y-axis
%     dx: dx of sampling points along x-axis
%     dx: dy of sampling points along y-axis
%     num_of_point_per_lay: [num_of_intfce] array, number of points of each layer
%     point_elev: {num_of_intfce} cell array, each elem is [nx,ny,num_of_intfce] array
%     
%     the meaning of possible vars ({num_of_intfce}[ny,nx, num_of_point_per_lay] array
%     attention: the order of array is [x,y, layer] as in md3grd file
%       val: only for one_component, keep target values
%       Vp
%       Vs
%       

%-- open file
fid = fopen(fnm_in,'r');

%-- 1st: media_type, value could be:
md.media_type = fscanf(fid, '%s', 1);

%-- 2nd: number of layer
md.num_of_intfce = fscanf(fid, '%d', 1);

%-- 3rd: number of points of each layer
md.num_of_point_per_lay = fscanf(fid, '%d', md.num_of_intfce);

%--- derived values
total_number_of_points_lay = sum(md.num_of_point_per_lay);

for n = 1 : md.num_of_intfce
  if n == 1
    start_index_of_layer(1) = 1;
  else
    start_index_of_layer(n) = sum(md.num_of_point_per_lay(1:n-1)) + 1;
  end
  end_index_of_layer(n)   = start_index_of_layer(n) + md.num_of_point_per_lay(n) - 1;
end

%-- 4rd
md.nx = fscanf(fid, '%d', 1);
md.ny = fscanf(fid, '%d', 1);
md.x0 = fscanf(fid, '%g', 1);
md.y0 = fscanf(fid, '%g', 1);
md.dx = fscanf(fid, '%g', 1);
md.dy = fscanf(fid, '%g', 1);

%-- rest
A = fscanf(fid, '%g');
%- attention: array order is [var, x,y,layer] in md3grd file

%-- create array
for n = 1 : md.num_of_intfce

    switch md.media_type

    %-- one component
    case 'one_component'
      if n == 1
        A = reshape(A,[2, md.nx, md.ny, total_number_of_points_lay]);
      end
      md.elev{n} = squeeze(A(1,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.val{n}  = squeeze(A(2,:,:, start_index_of_layer(n):end_index_of_layer(n)));

    %-- acoustic isotropic
    %   rho Vp
    case 'acoustic_isotropic'
      if n == 1
        A = reshape(A,[3, md.nx, md.ny, total_number_of_points_lay]);
      end
      md.elev{n}    = squeeze(A(1,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.density{n} = squeeze(A(2,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.Vp{n}      = squeeze(A(3,:,:, start_index_of_layer(n):end_index_of_layer(n)));

    %-- elastic isotropic
    %   rho Vp Vs
    case 'elastic_isotropic'
      if n == 1
        A = reshape(A,[4, md.nx, md.ny, total_number_of_points_lay]);
      end
      md.elev{n}    = squeeze(A(1,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.density{n} = squeeze(A(2,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.Vp{n}      = squeeze(A(3,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.Vs{n}      = squeeze(A(4,:,:, start_index_of_layer(n):end_index_of_layer(n)));

    %-- elastic vti, prem par
    %   rho Vph Vpv Vsh Vsv eta
    case 'elastic_vti_prem'
      if n == 1
        A = reshape(A,[7, md.nx, md.ny, total_number_of_points_lay]);
      end
      md.elev{n}    = squeeze(A(1,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.density{n} = squeeze(A(2,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.Vph{n}     = squeeze(A(3,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.Vpv{n}     = squeeze(A(4,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.Vsh{n}     = squeeze(A(5,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.Vsv{n}     = squeeze(A(6,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.eta{n}     = squeeze(A(7,:,:, start_index_of_layer(n):end_index_of_layer(n)));

    %-- elastic vti, thomsen par
    %   rho Vpv Vsv epsilon delta gamma
    case 'elastic_vti_thomsen'
      if n == 1
        A = reshape(A,[7, md.nx, md.ny, total_number_of_points_lay]);
      end
      md.elev{n}    = squeeze(A(1,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.density{n} = squeeze(A(2,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.Vpv{n}     = squeeze(A(3,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.Vsv{n}     = squeeze(A(4,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.epsilon{n} = squeeze(A(5,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.delta{n}   = squeeze(A(6,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.gamma{n}   = squeeze(A(7,:,:, start_index_of_layer(n):end_index_of_layer(n)));

    %-- elastic vti, thomsen par
    %   rho c11 c33 c55 c66 c13
    case 'elastic_vti_cij'
      if n == 1
        A = reshape(A,[7, md.nx, md.ny, total_number_of_points_lay]);
      end
      md.elev{n}    = squeeze(A(1,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.density{n} = squeeze(A(2,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C11{n}     = squeeze(A(3,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C33{n}     = squeeze(A(4,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C55{n}     = squeeze(A(5,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C66{n}     = squeeze(A(6,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C13{n}     = squeeze(A(7,:,:, start_index_of_layer(n):end_index_of_layer(n)));

    %-- elastic tti, thomsen par
    %   rho Vpv Vsv epsilon delta gamma azimuth dip
    case 'elastic_tti_thomsen'
      if n == 1
        A = reshape(A,[9, md.nx, md.ny, total_number_of_points_lay]);
      end
      md.elev{n}    = squeeze(A(1,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.density{n} = squeeze(A(2,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.Vpv{n}     = squeeze(A(3,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.Vsv{n}     = squeeze(A(4,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.epsilon{n} = squeeze(A(5,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.delta{n}   = squeeze(A(6,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.gamma{n}   = squeeze(A(7,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.azimuth{n} = squeeze(A(8,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.dip{n}     = squeeze(A(9,:,:, start_index_of_layer(n):end_index_of_layer(n)));

    %-- elastic tti, vti cij plus rotate
    %   rho c11 c33 c55 c66 c13 azimuth dip
    case 'elastic_tti_bond'
      if n == 1
        A = reshape(A,[9, md.nx, md.ny, total_number_of_points_lay]);
      end
      md.elev{n}    = squeeze(A(1,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.density{n} = squeeze(A(2,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C11{n}     = squeeze(A(3,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C33{n}     = squeeze(A(4,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C55{n}     = squeeze(A(5,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C66{n}     = squeeze(A(6,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C13{n}     = squeeze(A(7,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.azimuth{n} = squeeze(A(8,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.dip{n}     = squeeze(A(9,:,:, start_index_of_layer(n):end_index_of_layer(n)));

    %-- elastic aniso, cij
    %   rho c11 c12 c13 c14 c15 c16 c22 ...
    case 'elastic_aniso_cij'
      if n == 1
        A = reshape(A,[23, md.nx, md.ny, total_number_of_points_lay]);
      end
      md.elev{n}    = squeeze(A(1 ,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.density{n} = squeeze(A(2 ,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C11{n}     = squeeze(A(3 ,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C12{n}     = squeeze(A(4 ,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C13{n}     = squeeze(A(5 ,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C14{n}     = squeeze(A(6 ,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C15{n}     = squeeze(A(7 ,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C16{n}     = squeeze(A(8 ,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C22{n}     = squeeze(A(9 ,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C23{n}     = squeeze(A(10,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C24{n}     = squeeze(A(11,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C25{n}     = squeeze(A(12,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C26{n}     = squeeze(A(13,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C33{n}     = squeeze(A(14,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C34{n}     = squeeze(A(15,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C35{n}     = squeeze(A(16,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C36{n}     = squeeze(A(17,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C44{n}     = squeeze(A(18,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C45{n}     = squeeze(A(19,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C46{n}     = squeeze(A(20,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C55{n}     = squeeze(A(21,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C56{n}     = squeeze(A(22,:,:, start_index_of_layer(n):end_index_of_layer(n)));
      md.C66{n}     = squeeze(A(23,:,:, start_index_of_layer(n):end_index_of_layer(n)));
    end % swith media_type
  
end % n

fclose(fid);

end % function
