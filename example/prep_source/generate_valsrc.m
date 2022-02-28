%------------------------------------------------------------------------------
%-- Example of generating the value source file format 
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
flag_of_stf = 1;
% flag for source component and mechanism format
% 1st value: source components, 1(force), 2(momoment), 3(force+moment)
% 2nd value: mechanism format for moment source: 0 moment, 1 angle + mu + D + A
flag_of_cmp = 2;
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

dt_in = 0.02; 
nt_in = 51;   
% t_start
t_start = [1.0,0.0,2.0];
% t_start

fc = 2.;
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
% if mechism is mechanism_angle, in order input (strike dip rake u slip_rate area)
% moment_value
% % mechanism_angle
for i = 1 : nt_in
    strike(i) = 70.0;
    dip(i) = 30.0;
    rake(i) = 80.0;
    u(i) = 1.0e10;% u is shear modulus
    V(i) = 1.0;% unit is m/s
    A(i) = 1.0e6;%unit is meter^2
end
%==============================================================================
%-- write .valsrc file
%==============================================================================
anasrc_file = "test_source.valsrc";

fid=fopen(anasrc_file,'w'); % Output file name 
fprintf(fid,test_name); %test name
fprintf(fid,'%g\n',in_num_source);
fprintf(fid,'%g %g %g\n',flag_of_stf,dt_in,nt_in);
fprintf(fid,'%g %g\n',flag_of_cmp,mechanism_format);
fprintf(fid,'%g %g\n',flag_loc,flag_loc_z);
for i = 1 : in_num_source
    fprintf(fid,'%g %g %g\n',x_coord(i),y_coord(i),y_coord(i));
end

for i = 1 : in_num_source
  fprintf(fid,'%g\n',t_start(i));
  if(flag_of_cmp == 1 && mechanism_format == 0)
    for j = 1:nt_in
        fprintf(fid,'%g %g %g\n',x_vector(j),y_vector(j),z_vector(j));
    end
  end
  if(flag_of_cmp == 2 && mechanism_format == 0)
    for j = 1:nt_in
        fprintf(fid,'%g %g %g %g %g %g\n',Mxx(i),Myy(i),Mzz(i),Myz(i),Mxz(i),Mxy(i));
    end
  end
  if(flag_of_cmp == 2 && mechanism_format == 1)
    for j = 1:nt_in
        fprintf(fid,'%g %g %g %g %g %g\n',strike(i),dip(i),rake(i),u(i),V(i),A(i));
    end
  end
  %less used
  if(flag_of_cmp == 3)  
    for j = 1:nt_in
        fprintf(fid,'%g %g %g, %g %g %g %g %g %g\n',x_vector(j),y_vector(j),z_vector(j), Mxx(i),Myy(i),Mzz(i),Myz(i),Mxz(i),Mxy(i));
    end
  end
  
end
fclose(fid);

