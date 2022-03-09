function md = md3lay_import(fnm_in)

% This function exports structure model to md3lay file.
%  fnm_in: string, the md3lay file name;
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

%-- open file
fid = fopen(fnm_in,'r');

%-- 1st: media_type, value could be:
md.media_type = fscanf(fid, '%s', 1);

%-- 2nd: number of layer
md.num_of_intfce = fscanf(fid, '%d', 1);

%-- 3rd
md.nx = fscanf(fid, '%d', 1);
md.ny = fscanf(fid, '%d', 1);
md.x0 = fscanf(fid, '%g', 1);
md.y0 = fscanf(fid, '%g', 1);
md.dx = fscanf(fid, '%g', 1);
md.dy = fscanf(fid, '%g', 1);

%-- rest
A = fscanf(fid, '%g');
%- attention: array order is [3*var,x,y,layer] in md3lay file

%-- create array
for n = 1 : md.num_of_intfce

    switch md.media_type

    %-- one component
    case 'one_component'
      if n == 1
        A = reshape(A,[1+1*3, md.nx, md.ny, md.num_of_intfce]);
      end
      md.elev{n} = squeeze(A(1,:,:,n));
      md.val     {n} = squeeze(A(2,:,:,n));
      md.val_coef{n} = squeeze(A(3,:,:,n));
      md.val_pow {n} = squeeze(A(4,:,:,n));

    %-- acoustic isotropic
    %   rho Vp
    case 'acoustic_isotropic'
      if n == 1
        A = reshape(A,[1+2*3, md.nx, md.ny, md.num_of_intfce]);
      end
      md.elev{n} = squeeze(A(1,:,:,n));
      md.density     {n} = squeeze(A(2,:,:,n));
      md.density_coef{n} = squeeze(A(3,:,:,n));
      md.density_pow {n} = squeeze(A(4,:,:,n));
      md.Vp          {n} = squeeze(A(5,:,:,n));
      md.Vp_coef     {n} = squeeze(A(6,:,:,n));
      md.Vp_pow      {n} = squeeze(A(7,:,:,n));

    %-- elastic isotropic
    %   rho Vp Vs
    case 'elastic_isotropic'
      if n == 1
        A = reshape(A,[1+3*3, md.nx, md.ny, md.num_of_intfce]);
      end
      md.elev{n} = squeeze(A(1,:,:,n));
      md.density     {n} = squeeze(A(2,:,:,n));
      md.density_coef{n} = squeeze(A(3,:,:,n));
      md.density_pow {n} = squeeze(A(4,:,:,n));
      md.Vp          {n} = squeeze(A(5,:,:,n));
      md.Vp_coef     {n} = squeeze(A(6,:,:,n));
      md.Vp_pow      {n} = squeeze(A(7,:,:,n));
      md.Vs          {n} = squeeze(A(8,:,:,n));
      md.Vs_coef     {n} = squeeze(A(9,:,:,n));
      md.Vs_pow      {n} = squeeze(A(10,:,:,n));

    %-- elastic vti, prem par
    %   rho Vph Vpv Vsh Vsv eta
    case 'elastic_vti_prem'
      if n == 1
        A = reshape(A,[1+6*3, md.nx, md.ny, md.num_of_intfce]);
      end
      md.elev{n} = squeeze(A(1,:,:,n));
      md.density     {n} = squeeze(A( 2,:,:,n));
      md.density_coef{n} = squeeze(A( 3,:,:,n));
      md.density_pow {n} = squeeze(A( 4,:,:,n));
      md.Vph         {n} = squeeze(A( 5,:,:,n));
      md.Vph_coef    {n} = squeeze(A( 6,:,:,n));
      md.Vph_pow     {n} = squeeze(A( 7,:,:,n));
      md.Vpv         {n} = squeeze(A( 8,:,:,n));
      md.Vpv_coef    {n} = squeeze(A( 9,:,:,n));
      md.Vpv_pow     {n} = squeeze(A(10,:,:,n));
      md.Vsh         {n} = squeeze(A(11,:,:,n));
      md.Vsh_coef    {n} = squeeze(A(12,:,:,n));
      md.Vsh_pow     {n} = squeeze(A(13,:,:,n));
      md.Vsv         {n} = squeeze(A(14,:,:,n));
      md.Vsv_coef    {n} = squeeze(A(15,:,:,n));
      md.Vsv_pow     {n} = squeeze(A(16,:,:,n));
      md.eta         {n} = squeeze(A(17,:,:,n));
      md.eta_coef    {n} = squeeze(A(18,:,:,n));
      md.eta_pow     {n} = squeeze(A(19,:,:,n));

    %-- elastic vti, thomsen par
    %   rho Vpv Vsv epsilon delta gamma
    case 'elastic_vti_thomsen'
      if n == 1
        A = reshape(A,[1+6*3, md.nx, md.ny, md.num_of_intfce]);
      end
      md.elev{n} = squeeze(A(1,:,:,n));
      md.density     {n} = squeeze(A( 2,:,:,n));
      md.density_coef{n} = squeeze(A( 3,:,:,n));
      md.density_pow {n} = squeeze(A( 4,:,:,n));
      md.Vpv         {n} = squeeze(A( 5,:,:,n));
      md.Vpv_coef    {n} = squeeze(A( 6,:,:,n));
      md.Vpv_pow     {n} = squeeze(A( 7,:,:,n));
      md.Vsv         {n} = squeeze(A( 8,:,:,n));
      md.Vsv_coef    {n} = squeeze(A( 9,:,:,n));
      md.Vsv_pow     {n} = squeeze(A(10,:,:,n));
      md.epsilon     {n} = squeeze(A(11,:,:,n));
      md.epsilon_coef{n} = squeeze(A(12,:,:,n));
      md.epsilon_pow {n} = squeeze(A(13,:,:,n));
      md.delta       {n} = squeeze(A(14,:,:,n));
      md.delta_coef  {n} = squeeze(A(15,:,:,n));
      md.delta_pow   {n} = squeeze(A(16,:,:,n));
      md.gamma       {n} = squeeze(A(17,:,:,n));
      md.gamma_coef  {n} = squeeze(A(18,:,:,n));
      md.gamma_pow   {n} = squeeze(A(19,:,:,n));

    %-- elastic vti, thomsen par
    %   rho c11 c33 c55 c66 c13
    case 'elastic_vti_cij'
      if n == 1
        A = reshape(A,[1+6*3, md.nx, md.ny, md.num_of_intfce]);
      end
      md.elev{n} = squeeze(A(1,:,:,n));
      md.density     {n} = squeeze(A( 2,:,:,n));
      md.density_coef{n} = squeeze(A( 3,:,:,n));
      md.density_pow {n} = squeeze(A( 4,:,:,n));
      md.C11         {n} = squeeze(A( 5,:,:,n));
      md.C11_coef    {n} = squeeze(A( 6,:,:,n));
      md.C11_pow     {n} = squeeze(A( 7,:,:,n));
      md.C33         {n} = squeeze(A( 8,:,:,n));
      md.C33_coef    {n} = squeeze(A( 9,:,:,n));
      md.C33_pow     {n} = squeeze(A(10,:,:,n));
      md.C55         {n} = squeeze(A(11,:,:,n));
      md.C55_coef    {n} = squeeze(A(12,:,:,n));
      md.C55_pow     {n} = squeeze(A(13,:,:,n));
      md.C66         {n} = squeeze(A(14,:,:,n));
      md.C66_coef    {n} = squeeze(A(15,:,:,n));
      md.C66_pow     {n} = squeeze(A(16,:,:,n));
      md.C13         {n} = squeeze(A(17,:,:,n));
      md.C13_coef    {n} = squeeze(A(18,:,:,n));
      md.C13_pow     {n} = squeeze(A(19,:,:,n));

    %-- elastic tti, thomsen par
    %   rho Vpv Vsv epsilon delta gamma azimuth dip
    case 'elastic_tti_thomsen'
      if n == 1
        A = reshape(A,[1+8*3, md.nx, md.ny, md.num_of_intfce]);
      end
      md.elev{n} = squeeze(A(1,:,:,n));
      md.density     {n} = squeeze(A( 2,:,:,n));
      md.density_coef{n} = squeeze(A( 3,:,:,n));
      md.density_pow {n} = squeeze(A( 4,:,:,n));
      md.Vpv         {n} = squeeze(A( 5,:,:,n));
      md.Vpv_coef    {n} = squeeze(A( 6,:,:,n));
      md.Vpv_pow     {n} = squeeze(A( 7,:,:,n));
      md.Vsv         {n} = squeeze(A( 8,:,:,n));
      md.Vsv_coef    {n} = squeeze(A( 9,:,:,n));
      md.Vsv_pow     {n} = squeeze(A(10,:,:,n));
      md.epsilon     {n} = squeeze(A(11,:,:,n));
      md.epsilon_coef{n} = squeeze(A(12,:,:,n));
      md.epsilon_pow {n} = squeeze(A(13,:,:,n));
      md.delta       {n} = squeeze(A(14,:,:,n));
      md.delta_coef  {n} = squeeze(A(15,:,:,n));
      md.delta_pow   {n} = squeeze(A(16,:,:,n));
      md.gamma       {n} = squeeze(A(17,:,:,n));
      md.gamma_coef  {n} = squeeze(A(18,:,:,n));
      md.gamma_pow   {n} = squeeze(A(19,:,:,n));
      md.azimuth     {n} = squeeze(A(20,:,:,n));
      md.azimuth_coef{n} = squeeze(A(21,:,:,n));
      md.azimuth_pow {n} = squeeze(A(22,:,:,n));
      md.dip         {n} = squeeze(A(23,:,:,n));
      md.dip_coef    {n} = squeeze(A(24,:,:,n));
      md.dip_pow     {n} = squeeze(A(25,:,:,n));

    %-- elastic tti, vti cij plus rotate
    %   rho c11 c33 c55 c66 c13 azimuth dip
    case 'elastic_tti_bond'
      if n == 1
        A = reshape(A,[1+8*3, md.nx, md.ny, md.num_of_intfce]);
      end
      md.elev{n} = squeeze(A(1,:,:,n));
      md.density     {n} = squeeze(A( 2,:,:,n));
      md.density_coef{n} = squeeze(A( 3,:,:,n));
      md.density_pow {n} = squeeze(A( 4,:,:,n));
      md.C11         {n} = squeeze(A( 5,:,:,n));
      md.C11_coef    {n} = squeeze(A( 6,:,:,n));
      md.C11_pow     {n} = squeeze(A( 7,:,:,n));
      md.C33         {n} = squeeze(A( 8,:,:,n));
      md.C33_coef    {n} = squeeze(A( 9,:,:,n));
      md.C33_pow     {n} = squeeze(A(10,:,:,n));
      md.C55         {n} = squeeze(A(11,:,:,n));
      md.C55_coef    {n} = squeeze(A(12,:,:,n));
      md.C55_pow     {n} = squeeze(A(13,:,:,n));
      md.C66         {n} = squeeze(A(14,:,:,n));
      md.C66_coef    {n} = squeeze(A(15,:,:,n));
      md.C66_pow     {n} = squeeze(A(16,:,:,n));
      md.C13         {n} = squeeze(A(17,:,:,n));
      md.C13_coef    {n} = squeeze(A(18,:,:,n));
      md.C13_pow     {n} = squeeze(A(19,:,:,n));
      md.azimuth     {n} = squeeze(A(20,:,:,n));
      md.azimuth_coef{n} = squeeze(A(21,:,:,n));
      md.azimuth_pow {n} = squeeze(A(22,:,:,n));
      md.dip         {n} = squeeze(A(23,:,:,n));
      md.dip_coef    {n} = squeeze(A(24,:,:,n));
      md.dip_pow     {n} = squeeze(A(25,:,:,n));

    %-- elastic aniso, cij
    %   rho c11 c12 c13 c14 c15 c16 c22 ...
    case 'elastic_aniso_cij'
      if n == 1
        A = reshape(A,[1+22*3, md.nx, md.ny, md.num_of_intfce]);
      end
      md.elev{n} = squeeze(A(1,:,:,n));
      md.density     {n} = squeeze(A( 2,:,:,n));
      md.density_coef{n} = squeeze(A( 3,:,:,n));
      md.density_pow {n} = squeeze(A( 4,:,:,n));
      md.C11         {n} = squeeze(A( 5,:,:,n));
      md.C11_coef    {n} = squeeze(A( 6,:,:,n));
      md.C11_pow     {n} = squeeze(A( 7,:,:,n));
      md.C12         {n} = squeeze(A( 8,:,:,n));
      md.C12_coef    {n} = squeeze(A( 9,:,:,n));
      md.C12_pow     {n} = squeeze(A(10,:,:,n));
      md.C13         {n} = squeeze(A(11,:,:,n));
      md.C13_coef    {n} = squeeze(A(12,:,:,n));
      md.C13_pow     {n} = squeeze(A(13,:,:,n));
      md.C14         {n} = squeeze(A(14,:,:,n));
      md.C14_coef    {n} = squeeze(A(15,:,:,n));
      md.C14_pow     {n} = squeeze(A(16,:,:,n));
      md.C15         {n} = squeeze(A(17,:,:,n));
      md.C15_coef    {n} = squeeze(A(18,:,:,n));
      md.C15_pow     {n} = squeeze(A(19,:,:,n));
      md.C16         {n} = squeeze(A(20,:,:,n));
      md.C16_coef    {n} = squeeze(A(21,:,:,n));
      md.C16_pow     {n} = squeeze(A(22,:,:,n));
      md.C22         {n} = squeeze(A(23,:,:,n));
      md.C22_coef    {n} = squeeze(A(25,:,:,n));
      md.C22_pow     {n} = squeeze(A(25,:,:,n));
      md.C23         {n} = squeeze(A(26,:,:,n));
      md.C23_coef    {n} = squeeze(A(27,:,:,n));
      md.C23_pow     {n} = squeeze(A(28,:,:,n));
      md.C24         {n} = squeeze(A(29,:,:,n));
      md.C24_coef    {n} = squeeze(A(30,:,:,n));
      md.C24_pow     {n} = squeeze(A(31,:,:,n));
      md.C25         {n} = squeeze(A(32,:,:,n));
      md.C25_coef    {n} = squeeze(A(33,:,:,n));
      md.C25_pow     {n} = squeeze(A(34,:,:,n));
      md.C26         {n} = squeeze(A(35,:,:,n));
      md.C26_coef    {n} = squeeze(A(36,:,:,n));
      md.C26_pow     {n} = squeeze(A(37,:,:,n));
      md.C33         {n} = squeeze(A(38,:,:,n));
      md.C33_coef    {n} = squeeze(A(39,:,:,n));
      md.C33_pow     {n} = squeeze(A(40,:,:,n));
      md.C34         {n} = squeeze(A(41,:,:,n));
      md.C34_coef    {n} = squeeze(A(42,:,:,n));
      md.C34_pow     {n} = squeeze(A(43,:,:,n));
      md.C35         {n} = squeeze(A(44,:,:,n));
      md.C35_coef    {n} = squeeze(A(45,:,:,n));
      md.C35_pow     {n} = squeeze(A(46,:,:,n));
      md.C36         {n} = squeeze(A(47,:,:,n));
      md.C36_coef    {n} = squeeze(A(48,:,:,n));
      md.C36_pow     {n} = squeeze(A(49,:,:,n));
      md.C44         {n} = squeeze(A(50,:,:,n));
      md.C44_coef    {n} = squeeze(A(51,:,:,n));
      md.C44_pow     {n} = squeeze(A(52,:,:,n));
      md.C45         {n} = squeeze(A(53,:,:,n));
      md.C45_coef    {n} = squeeze(A(54,:,:,n));
      md.C45_pow     {n} = squeeze(A(55,:,:,n));
      md.C46         {n} = squeeze(A(56,:,:,n));
      md.C46_coef    {n} = squeeze(A(57,:,:,n));
      md.C46_pow     {n} = squeeze(A(58,:,:,n));
      md.C55         {n} = squeeze(A(59,:,:,n));
      md.C55_coef    {n} = squeeze(A(60,:,:,n));
      md.C55_pow     {n} = squeeze(A(61,:,:,n));
      md.C56         {n} = squeeze(A(62,:,:,n));
      md.C56_coef    {n} = squeeze(A(63,:,:,n));
      md.C56_pow     {n} = squeeze(A(64,:,:,n));
      md.C66         {n} = squeeze(A(65,:,:,n));
      md.C66_coef    {n} = squeeze(A(66,:,:,n));
      md.C66_pow     {n} = squeeze(A(67,:,:,n));
    end % swith media_type

end % n

fclose(fid);

end % function
