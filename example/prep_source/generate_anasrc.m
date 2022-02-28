%------------------------------------------------------------------------------
%-- Example of generating the analytic source file format
%------------------------------------------------------------------------------
% NOTE Only matlab 2016 and above versions are supported
clear all;
close all;
clc;
test_name = "test_event_1\n"; %test name
in_num_source = 3;
%flag for stf
%1st value : 0 analytic stf or 1 discrete values
%for analytical, 2nd value is time length
%for discrete,  2nd value is dt and 3rd is nt
flag_of_stf = 0;
t_len = 4.0; %source time window length
% flag for source component and mechanism format
% 1st value: source components, 1(force), 2(momoment), 3(force+moment)
% 2nd value: mechanism format for moment source: 0 moment, 1 angle + mu + D + A
flag_of_cmp = 0;
mechanism_format = 0; 
% flag for location
% 1st value: 0 computational coordiate(index), 1 physical coordinate(coords)
% 2nd value: third coordinate is 0 axis or 1 depth
flag_loc = 1;
flag_loc_z = 0;
% coords
x_coord = [1000.0, 2000.0, 3000.0 ];
y_coord = [2200.0, 3050.0, 5000.0 ];
z_coord = [3200.0, 2300.0, 1000.0 ];


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


% source para fc and t0
fc = [2.0, 2.0, 2.0];
t0 = [0.5, 0.5, 0.5];

% t_start
t_start = [1.0,0.0,2.0];
%sorce type name
%force source include type  (ricker, gaussian)
%moment source include type  (ricker_deriv, gaussian_deriv)
wavelet_name = ["ricker","ricker","ricker"];

%==============================================================================
%-- write .anasrc file
%==============================================================================
anasrc_file = "test_source.anasrc";

fid=fopen(anasrc_file,'w'); % Output file name 
fprintf(fid,test_name); 
fprintf(fid,'%g\n',in_num_source);
fprintf(fid,'%g %g\n',flag_of_stf,t_len);
fprintf(fid,'%g %g\n',flag_of_cmp,mechanism_format);
fprintf(fid,'%g %g\n',flag_loc,flag_loc_z);
for i = 1 : in_num_source
    fprintf(fid,'%g %g %g\n',x_coord(i),y_coord(i),y_coord(i));
end
for i = 1 : in_num_source
  fprintf(fid,'%g %s %g %g\n',t_start(i),wavelet_name(i),fc(i),t0(i));
  if(flag_of_cmp == 0)
      fprintf(fid,'%g %g %g\n',x_vector(i),y_vector(i),y_vector(i));
  end
  if(flag_of_cmp == 1 && mechanism_format == 0)
      fprintf(fid,'%g %g %g %g %g %g\n',Mxx(i),Myy(i),Mzz(i),Myz(i),Mxz(i),Mxy(i));
  end
  if(flag_of_cmp == 1 && mechanism_format == 1)
      fprintf(fid,'%g %g %g %g %g %g\n',strike(i),dip(i),rake(i),u(i),D(i),A(i));
  end
end

fclose(fid);
