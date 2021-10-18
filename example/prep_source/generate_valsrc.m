%------------------------------------------------------------------------------
%-- Example of generating the value source file format 
%------------------------------------------------------------------------------
% NOTE Only matlab 2016 and above versions are supported
clear all;
close all;
clc;
test_name = "test_event_1\n"; %test name
num_force_src = 2;
num_moment_src = 1;
dt_in = 0.02; %
nt_in = 51;   %
% coords
x_force = [1000.0, 2000.0];
y_force = [2200.0, 3050.0];
z_force = [3200.0, 2300.0];
x_moment = [2000.0];
y_moment = [3200.0];
z_moment = [4200.0];

%include two type (moment_tensor, mechanism_angle)
moment_wavelet_mechism = ["moment_tensor\n"];
% moment_wavelet_mechism = ["mechanism_angle\n"];

% t_start
t_force_start = [1.0,0.0];
t_moment_start = [1.0];
fc = 2.0;
t0 = 0.5;
% force_value 
% Note component not need equal. here equal is just for easy test
for i = 1 : nt_in
    t = (i-1) * dt_in;
    x_vector(i) = stf_rickerderiv(t,fc,t0);
    y_vector(i) = stf_rickerderiv(t,fc,t0);
    z_vector(i) = stf_rickerderiv(t,fc,t0);
    Mxx(i) = stf_rickerderiv(t,fc,t0);
    Myy(i) = stf_rickerderiv(t,fc,t0);
    Mzz(i) = stf_rickerderiv(t,fc,t0);
    Myz(i) = stf_rickerderiv(t,fc,t0);
    Mxz(i) = stf_rickerderiv(t,fc,t0);
    Mxy(i) = stf_rickerderiv(t,fc,t0);
end
% if mechism is moment_tensor, in oeder input (Mxx, Myy, Mzz, Myz, Mxz, Mxy).
% if mechism is mechanism_angle, in order input (strike dip rake u slip area)
% moment_value
% % mechanism_angle
for i = 1 : nt_in
    strike(i) = 70.0;
    dip(i) = 30.0;
    rake(i) = 80.0;
    u(i) = 1.0e10;% u is shear modulus
    D(i) = 1.0;%unit is meter
    A(i) = 1.0e6;%unit is meter^2
end
%==============================================================================
%-- write .valsrc file
%==============================================================================
anasrc_file = "test_source.valsrc";

fid=fopen(anasrc_file,'w'); % Output file name 
fprintf(fid,test_name); %test name
fprintf(fid,'%g %g\n',num_force_src,num_moment_src);
fprintf(fid,'%g %g\n',dt_in,nt_in);  %source time window length

for i = 1 : num_force_src
    fprintf(fid,'%g %g %g\n',x_force(i),y_force(i),z_force(i));
end
for i = 1 : num_moment_src
  fprintf(fid,'%g %g %g\n',x_moment(i),y_moment(i),z_moment(i));
end

for i = 1 : num_force_src
  fprintf(fid,'%g\n',t_force_start(i));
  for j = 1:nt_in
      fprintf(fid,'%g %g %g\n',x_vector(j),y_vector(j),z_vector(j));
  end
end

for i = 1 : num_moment_src
  fprintf(fid, '%g\n',t_moment_start(i));
  fprintf(fid,moment_wavelet_mechism(i));
  if (moment_wavelet_mechism(i) == 'moment_tensor\n')
  for j = 1:nt_in
      fprintf(fid,'%g %g %g %g %g %g\n',Mxx(j),Myy(j),Mzz(j), ...
          Myz(j),Mxz(j),Mxy(j));
  end
  end
  if (moment_wavelet_mechism(i) == 'mechanism_angle\n')
  for j = 1:nt_in
      fprintf(fid,'%g %g %g %g %g %g\n',strike(j),dip(j),rake(j), ...
          u(j),D(j),A(j));
  end
  end
end
fclose(fid);

