function md3lay_export(fnm_ou, md)

% This function exports structure model to md3lay file.
%  fnm_ou: string, the output file name;
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
%     elev: {num_of_intfce} cell array, each elem is [nx,ny,num_of_intfce] array
%     
%     the meaning of possible vars ({num_of_intfce}[ny,nx, num_of_point_per_lay] array
%       val: only for one_component, keep target values
%       val_coef: coefficient of polynomail
%       val_pow: exponent of polynomail in relative depth
%       Vp
%       Vs
%       

%-- create output file
fid = fopen(fnm_ou,'w');

%-- 1st: media_type, value could be:
fprintf(fid, '%s\n',md.media_type);

%-- 2nd: number of layer
fprintf(fid, '%d\n', md.num_of_intfce);

%-- 3rd
fprintf(fid, '%d %d %f %f %f %f\n', md.nx, md.ny, md.x0, md.y0, md.dx, md.dy);

%-- rest
  for ilay = 1 : md.num_of_intfce
      disp([num2str(ilay), 'th-layer of total ', num2str(md.num_of_intfce), ' layers'])
      for j = 1 : md.ny
          for i = 1 : md.nx
              % elevation
            	fprintf(fid, '%g',  md.elev{ilay}(i,j));

              switch md.media_type

              %-- one component
              case 'one_component'
            	  fprintf(fid, ' %g %g %g', ...
                              md.val     {ilay}(i,j), ...
                              md.val_coef{ilay}(i,j), ...
                              md.val_pow {ilay}(i,j)  ...
                              );

              %-- acoustic isotropic
              %   rho Vp
              case 'acoustic_isotropic'
            	  fprintf(fid, ' %g %g %g', ...
                              md.density     {ilay}(i,j), ...
                              md.density_coef{ilay}(i,j), ...
                              md.density_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.Vp     {ilay}(i,j), ...
                              md.Vp_coef{ilay}(i,j), ...
                              md.Vp_pow {ilay}(i,j) );

              %-- elastic isotropic
              %   rho Vp Vs
              case 'elastic_isotropic'
            	  fprintf(fid, ' %g %g %g', ...
                              md.density     {ilay}(i,j), ...
                              md.density_coef{ilay}(i,j), ...
                              md.density_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.Vp     {ilay}(i,j), ...
                              md.Vp_coef{ilay}(i,j), ...
                              md.Vp_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.Vs     {ilay}(i,j), ...
                              md.Vs_coef{ilay}(i,j), ...
                              md.Vs_pow {ilay}(i,j) );

              %-- elastic vti, prem par
              %   rho Vph Vpv Vsh Vsv eta
              case 'elastic_vti_prem'
            	  fprintf(fid, ' %g %g %g', ...
                              md.density     {ilay}(i,j), ...
                              md.density_coef{ilay}(i,j), ...
                              md.density_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.Vph     {ilay}(i,j), ...
                              md.Vph_coef{ilay}(i,j), ...
                              md.Vph_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.Vpv     {ilay}(i,j), ...
                              md.Vpv_coef{ilay}(i,j), ...
                              md.Vpv_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.Vsh     {ilay}(i,j), ...
                              md.Vsh_coef{ilay}(i,j), ...
                              md.Vsh_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.Vsv     {ilay}(i,j), ...
                              md.Vsv_coef{ilay}(i,j), ...
                              md.Vsv_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.eta     {ilay}(i,j), ...
                              md.eta_coef{ilay}(i,j), ...
                              md.eta_pow {ilay}(i,j) );

              %-- elastic vti, thomsen par
              %   rho Vpv Vsv epsilon delta gamma
              case 'elastic_vti_thomsen'
            	  fprintf(fid, ' %g %g %g', ...
                              md.density     {ilay}(i,j), ...
                              md.density_coef{ilay}(i,j), ...
                              md.density_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.Vpv     {ilay}(i,j), ...
                              md.Vpv_coef{ilay}(i,j), ...
                              md.Vpv_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.Vsv     {ilay}(i,j), ...
                              md.Vsv_coef{ilay}(i,j), ...
                              md.Vsv_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.epsilon     {ilay}(i,j), ...
                              md.epsilon_coef{ilay}(i,j), ...
                              md.epsilon_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.delta     {ilay}(i,j), ...
                              md.delta_coef{ilay}(i,j), ...
                              md.delta_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.gamma     {ilay}(i,j), ...
                              md.gamma_coef{ilay}(i,j), ...
                              md.gamma_pow {ilay}(i,j) );

              %-- elastic vti, thomsen par
              %   rho c11 c33 c55 c66 c13
              case 'elastic_vti_cij'
            	  fprintf(fid, ' %g %g %g', ...
                              md.density     {ilay}(i,j), ...
                              md.density_coef{ilay}(i,j), ...
                              md.density_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C11     {ilay}(i,j), ...
                              md.C11_coef{ilay}(i,j), ...
                              md.C11_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C33     {ilay}(i,j), ...
                              md.C33_coef{ilay}(i,j), ...
                              md.C33_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C55     {ilay}(i,j), ...
                              md.C55_coef{ilay}(i,j), ...
                              md.C55_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C66     {ilay}(i,j), ...
                              md.C66_coef{ilay}(i,j), ...
                              md.C66_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C13     {ilay}(i,j), ...
                              md.C13_coef{ilay}(i,j), ...
                              md.C13_pow {ilay}(i,j) );

              %-- elastic tti, thomsen par
              %   rho Vpv Vsv epsilon delta gamma azimuth dip
              case 'elastic_tti_thomsen'
            	  fprintf(fid, ' %g %g %g', ...
                              md.density     {ilay}(i,j), ...
                              md.density_coef{ilay}(i,j), ...
                              md.density_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.Vpv     {ilay}(i,j), ...
                              md.Vpv_coef{ilay}(i,j), ...
                              md.Vpv_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.Vsv     {ilay}(i,j), ...
                              md.Vsv_coef{ilay}(i,j), ...
                              md.Vsv_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.epsilon     {ilay}(i,j), ...
                              md.epsilon_coef{ilay}(i,j), ...
                              md.epsilon_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.delta     {ilay}(i,j), ...
                              md.delta_coef{ilay}(i,j), ...
                              md.delta_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.gamma     {ilay}(i,j), ...
                              md.gamma_coef{ilay}(i,j), ...
                              md.gamma_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.azimuth     {ilay}(i,j), ...
                              md.azimuth_coef{ilay}(i,j), ...
                              md.azimuth_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.dip     {ilay}(i,j), ...
                              md.dip_coef{ilay}(i,j), ...
                              md.dip_pow {ilay}(i,j) );

              %-- elastic tti, vti cij plus rotate
              %   rho c11 c33 c55 c66 c13 azimuth dip
              case 'elastic_tti_bond'
            	  fprintf(fid, ' %g %g %g', ...
                              md.density     {ilay}(i,j), ...
                              md.density_coef{ilay}(i,j), ...
                              md.density_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C11     {ilay}(i,j), ...
                              md.C11_coef{ilay}(i,j), ...
                              md.C11_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C33     {ilay}(i,j), ...
                              md.C33_coef{ilay}(i,j), ...
                              md.C33_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C55     {ilay}(i,j), ...
                              md.C55_coef{ilay}(i,j), ...
                              md.C55_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C66     {ilay}(i,j), ...
                              md.C66_coef{ilay}(i,j), ...
                              md.C66_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C13     {ilay}(i,j), ...
                              md.C13_coef{ilay}(i,j), ...
                              md.C13_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.azimuth     {ilay}(i,j), ...
                              md.azimuth_coef{ilay}(i,j), ...
                              md.azimuth_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.dip     {ilay}(i,j), ...
                              md.dip_coef{ilay}(i,j), ...
                              md.dip_pow {ilay}(i,j) );

              %-- elastic aniso, cij
              %   rho c11 c12 c13 c14 c15 c16 c22 ...
              case 'elastic_aniso_cij'
            	  fprintf(fid, ' %g %g %g', ...
                              md.density     {ilay}(i,j), ...
                              md.density_coef{ilay}(i,j), ...
                              md.density_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C11     {ilay}(i,j), ...
                              md.C11_coef{ilay}(i,j), ...
                              md.C11_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C12     {ilay}(i,j), ...
                              md.C12_coef{ilay}(i,j), ...
                              md.C12_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C13     {ilay}(i,j), ...
                              md.C13_coef{ilay}(i,j), ...
                              md.C13_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C14     {ilay}(i,j), ...
                              md.C14_coef{ilay}(i,j), ...
                              md.C14_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C15     {ilay}(i,j), ...
                              md.C15_coef{ilay}(i,j), ...
                              md.C15_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C16     {ilay}(i,j), ...
                              md.C16_coef{ilay}(i,j), ...
                              md.C16_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C22     {ilay}(i,j), ...
                              md.C22_coef{ilay}(i,j), ...
                              md.C22_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C23     {ilay}(i,j), ...
                              md.C23_coef{ilay}(i,j), ...
                              md.C23_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C24     {ilay}(i,j), ...
                              md.C24_coef{ilay}(i,j), ...
                              md.C24_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C25     {ilay}(i,j), ...
                              md.C25_coef{ilay}(i,j), ...
                              md.C25_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C26     {ilay}(i,j), ...
                              md.C26_coef{ilay}(i,j), ...
                              md.C26_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C33     {ilay}(i,j), ...
                              md.C33_coef{ilay}(i,j), ...
                              md.C33_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C34     {ilay}(i,j), ...
                              md.C34_coef{ilay}(i,j), ...
                              md.C34_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C35     {ilay}(i,j), ...
                              md.C35_coef{ilay}(i,j), ...
                              md.C35_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C36     {ilay}(i,j), ...
                              md.C36_coef{ilay}(i,j), ...
                              md.C36_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C44     {ilay}(i,j), ...
                              md.C44_coef{ilay}(i,j), ...
                              md.C44_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C45     {ilay}(i,j), ...
                              md.C45_coef{ilay}(i,j), ...
                              md.C45_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C46     {ilay}(i,j), ...
                              md.C46_coef{ilay}(i,j), ...
                              md.C46_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C55     {ilay}(i,j), ...
                              md.C55_coef{ilay}(i,j), ...
                              md.C55_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C56     {ilay}(i,j), ...
                              md.C56_coef{ilay}(i,j), ...
                              md.C56_pow {ilay}(i,j) );
            	  fprintf(fid, ' %g %g %g', ...
                              md.C66     {ilay}(i,j), ...
                              md.C66_coef{ilay}(i,j), ...
                              md.C66_pow {ilay}(i,j) );
              end % swith media_type
              
              % return
            	fprintf(fid, '\n');
          end % i
      end % j
  end % ilay

fclose(fid);

end % function
