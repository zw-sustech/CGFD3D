%------------------------------------------------------------------------------
%-- Example of generating the value source file format 
%------------------------------------------------------------------------------

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

% force_value 
x_vector = [4294882099200.00;8542101700608.00;16393744416768.0;30346210967552.0;
           54152927379456.0;93104052895744.0;154109701259264;245370726645760;
           375383295787008;551050343874560;774830152482816;1.04114704875520e+15;
           1.33270258214502e+15;1.61778586825523e+15;1.84997609524429e+15;1.97156474034586e+15;
           1.92146381262029e+15;1.64720464940237e+15;1.11916741138842e+15;343776363020288;
           -628307611615232;-1.70257469669376e+15;-2.75141832684339e+15;-3.63437206641050e+15;
           -4.22395564890522e+15;-4.43113493915238e+15;-4.22395618577613e+15;-3.63437448232960e+15;
           -2.75141832684339e+15;-1.70257644152422e+15;-628311101276160;343776363020288;
           1.11916653897318e+15;1.64720384409600e+15;1.92146381262029e+15;1.97156500878131e+15;
           1.84997649789747e+15;1.61778586825523e+15;1.33270258214502e+15;1.04114738429952e+15;
           774830957789184;551050343874560;375383698440192;245371162853376;
           154109701259264;93104153559040.0;54152990294016.0;30346210967552.0;
           16393771679744.0;8542116380672.00;4294882099200.00];
% Note component not need equal. here equal is just for easy test
y_vector = x_vector;
z_vector = x_vector;
% if mechism is moment_tensor, in oeder input (Mxx, Myy, Mzz, Myz, Mxz, Mxy).
% if mechism is mechanism_angle, in order input (strike dip rake u slip area)
% moment_value
Mxx = x_vector;
Myy = x_vector;
Mzz = x_vector;
Myz = x_vector;
Mxz = x_vector;
Mxy = x_vector;
% % mechanism_angle
% strike= x_vector;
% dip = x_vector;
% rake = x_vector;
% u = x_vector;% u is shear modulus
% D = x_vector;%unit is meter
% A = x_vector;%unit is meter^2
%==============================================================================
%-- write .valsrc file
%==============================================================================
anasrc_file = 'test_source.valsrc';

fid=fopen(anasrc_file,'w'); % Output file name 
fprintf(fid,test_name); %test name
fprintf(fid,'%d %d\n',num_force_src,num_moment_src);
fprintf(fid,'%.5f %d\n',dt_in,nt_in);  %source time window length

for i = 1 : num_force_src
    fprintf(fid,'%.5f %.5f %.5f\n',x_force(i),y_force(i),z_force(i));
end
for i = 1 : num_moment_src
  fprintf(fid,'%.5f %.5f %.5f\n',x_moment(i),y_moment(i),z_moment(i));
end

for i = 1 : num_force_src
  fprintf(fid,'%.5f\n',t_force_start(i));
  for j = 1:nt_in
      fprintf(fid,"%.5f %.5f %.5f\n",x_vector(j),y_vector(j),z_vector(j));
  end
end

for i = 1 : num_moment_src
  fprintf(fid, '%.5f\n',t_moment_start(i));
  fprintf(fid,moment_wavelet_mechism(i));
  if (moment_wavelet_mechism(i) == "moment_tensor\n")
  for j = 1:nt_in
      fprintf(fid,"%.5f %.5f %.5f %.5f, %.5f %.5f\n",Mxx(j),Myy(j),Mzz(j), ...
          Myz(j),Mxz(j),Mxy(j));
  end
  end
  if (moment_wavelet_mechism(i) == "mechanism_angle\n")
  for j = 1:nt_in
      fprintf(fid,"%.5f %.5f %.5f %.5f, %.5f %.5f\n",strike(j),dip(j),rake(j), ...
          u(j),D(j),A(j));
  end
  end
end


fclose(fid);
