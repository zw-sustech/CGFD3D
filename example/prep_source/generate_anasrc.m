%------------------------------------------------------------------------------
%-- Example of generating the analytic source file format
%------------------------------------------------------------------------------
% NOTE Only matlab 2016 and above versions are supported
clear all;
close all;
clc;
test_name = "test_event_1\n"; %test name
num_force_src = 3;
num_moment_src = 3;
t_len = 4.0; %source time window length
% coords
x_force = [1000.0, 2000.0, 3000.0 ];
y_force = [2200.0, 3050.0, 5000.0 ];
z_force = [3200.0, 2300.0, 1000.0 ];
x_moment = [2000.0, 1020.0, 2000.0 ];
y_moment = [3200.0, 2050.0, 3000.0 ];
z_moment = [4200.0, 3300.0, 2300.0 ];
% force_vector and moment_tensor
x_vector = [1e16, 1e16, 1e16];
y_vector = [1e16, 1e16, 1e16];
z_vector = [1e15, 1e15, 1e15];
% if mechism is moment_tensor, in oeder input (Mxx, Myy, Mzz, Myz, Mxz, Mxy).
% if mechism is mechanism_angle, in order input (strike dip rake u slip area)
Mxx = [1e15, 1e15, 1e15];
Myy = [1e15, 1e15, 1e15];
Mzz = [1e15, 1e15, 1e15];
Myz = [1e15, 1e15, 1e15];
Mxz = [1e15, 1e15, 1e15];
Mxy = [1e15, 1e15, 1e15];

strike = [80.0, 70.0, 50.0];
dip = [70.0, 60.0, 50.0];
rake = [30.0, 40.0, 50.0];
u = [1e10, 1e10, 1e10]; % u is shear modulus
D = [1.0, 1.1, 1.2]; %unit is meter
A = [1.0e8, 1.0e8, 1.0e8]; %unit is meter^2

%include two type (moment_tensor, mechanism_angle)
moment_wavelet_mechism = ["moment_tensor\n","moment_tensor\n","moment_tensor\n"];
% moment_wavelet_mechism = ["mechanism_angle\n","mechanism_angle\n","mechanism_angle\n"];

% source para fc and t0
fc_force = [2.0, 2.0, 2.0];
t0_force = [0.5, 0.5, 0.5];
fc_moment = [2.0, 2.0, 2.0];
t0_moment = [0.5, 0.5, 0.5];
% t_start
t_force_start = [1.0,0.0,2.0];
t_moment_start = [0.0,1.0,2.0];
%sorce type name
%force source include type  (ricker, gaussian)
force_wavelet_name = ["ricker\n","ricker\n","ricker\n"];
%moment source include type  (ricker_deriv, gaussian_deriv)
moment_wavelet_name = ["ricker_deriv\n","ricker_deriv\n","ricker_deriv\n"]; 

%==============================================================================
%-- write .anasrc file
%==============================================================================
anasrc_file = "test_source.anasrc";

fid=fopen(anasrc_file,'w'); % Output file name 
fprintf(fid,test_name); %test name
fprintf(fid,'%g %g\n',num_force_src,num_moment_src);
fprintf(fid,'%g\n',t_len);  %source time window length
for i = 1 : num_force_src
  fprintf(fid,'%g  %g  %g\n',x_force(i),y_force(i),z_force(i));
  fprintf(fid,'%g  %g  %g\n',x_vector(i),y_vector(i),z_vector(i));
  fprintf(fid, force_wavelet_name(i));
  fprintf(fid,'%g %g\n',fc_force(i),t0_force(i));
  fprintf(fid,'%g\n',t_force_start(i));
end
for i = 1 : num_moment_src
  fprintf(fid,'%g  %g  %g\n',x_moment(i),y_moment(i),z_moment(i));
  fprintf(fid,moment_wavelet_mechism(i));
  if (moment_wavelet_mechism(i) == "moment_tensor\n")
  fprintf(fid,'%g %g %g %g %g %g\n',Mxx(i),Myy(i),Mzz(i),Myz(i),Mxz(i),Mxy(i));
  end
  if (moment_wavelet_mechism(i) == "mechanism_angle\n")
  fprintf(fid,'%g %g %g %g %g %g\n',strike(i),dip(i),rake(i),u(i),D(i),A(i));
  end
  fprintf(fid,moment_wavelet_name(i));
  fprintf(fid, '%g %g\n',fc_moment(i),t0_moment(i));
  fprintf(fid, '%g\n',t_moment_start(i));
end
fclose(fid);
