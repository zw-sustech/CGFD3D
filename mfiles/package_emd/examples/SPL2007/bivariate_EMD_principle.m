%bivariate_EMD_principle.m 
%shows principle of the bivariate EMD extension
%reproduces Fig. 1 in "Bivariate Empirical Mode Decomposition", G. Rilling,
%P. Flandrin, P. Goncalves and J. M. Lilly, IEEE Signal Processing Letters
%
%G. Rilling 3/2007 email:  gabriel.rilling@ens-lyon.fr
function bivariate_EMD_principle
N = 8192; % number of points
tmax = 40;
Ndirs1 = 2; % used to display the envelope curves
Ndirs2 = 8; % used to render a nice 3D tube envelope

% display parameters
A = [0,tmax,-10,10,-10,10];
L = 2; % thickness of the envelope curves
gray = [.5,.5,.5];
target_faces_number = 1e3; % ~number of faces of the tube patch object. reduce/increase this number to get faster/nicer rendering
propNames = {
    'AmbientStrength'
    'DiffuseStrength'
    'SpecularStrength'
    'SpecularExponent'
    'SpecularColorReflectance'
    'FaceColor'
    'FaceAlpha'
    'EdgeColor'
  };
propVals = {0,0,3,7,1,gray,.3,'none'};
Axprops = {
  'XTickLabel'
  'YTickLabel'
  'ZTickLabel'
  'CameraPosition'
  'CameraViewAngle'
  };
CameraPosition = [-69.1664 -166.3868 17.2944];
CameraViewAngle = 7.4163;
Axpropvals = {[],[],[],CameraPosition,CameraViewAngle};
label_pos = [.07,.07,.09];
subplot_label = @(label)text(A(1)+(A(2)-A(1))*label_pos(1),A(3)+(A(4)-A(3))*(1-label_pos(2)),A(5)+(A(6)-A(5))*(1-label_pos(3)),['(',label,')']);

% abscissa vector
absc = linspace(0,tmax,N);
dt = mean(diff(absc));
% extended abscissa to avoid edge effects on the final plot
absc_calc = [-5:dt:-dt,absc,tmax+dt:dt:tmax+5];
tt = find(absc_calc>=0 & absc_calc <=tmax);

% fast rotating component
fm1 = @(t)1+.4*sin(2*pi/tmax*t);
am1 = @(t)1+.5*cos(pi/(2*tmax)*t);
dirstretch1 = @(t) 3*exp(i*t*pi/(1.5*tmax));
x1 = dirstretch(am1(absc_calc).*fmodany(dt*fm1(absc_calc).').',dirstretch1(absc_calc));

% slow rotating component
fm2 = @(t) 0.0382+t-t;
am2 = @(t) .2*(1+t);
x2 = am2(absc_calc).*fmodany(dt*fm2(absc_calc).').'; 

% composite signal
x = x1+x2;

% compute the envelope curves
% for a large number of directions to get a smooth tube envelope

[env2] = cenvelope(x,Ndirs2);
% env2 = flipud([emax;emin]);
[env1,moy] = cenvelope(x,Ndirs1);
% env1 = flipud([emax;emin]);


% make MATLAB render an extended tube envelope in order to get
% good normals at the boundaries of the final patch object
[m,n] = size(env2);
env = [env2;env2(1:3,:)];
step = max(round(numel(env)/target_faces_number),1);
inds = [fliplr(tt(1)-step:-step:1),tt(1):step:n];
if inds(end)<n
    inds = [inds,n];
end
inds2 = find(inds>=tt(1) & inds<= tt(end));
PYZ = env(:,inds);
PX = repmat(absc_calc(inds),m+3,1);
tmp = figure;
h = surf(PX,real(PYZ),imag(PYZ));
VN = get(h,'VertexNormals');
VN = VN(2:end-1,inds2,:);
PYZ = PYZ(2:end-1,inds2);
PX = PX(2:end-1,inds2);
delete(h);
close(tmp);


figure('defaultAxesNextplot','add','renderer','opengl','name','Principle of the Bivariate EMD');
% NOTE: openGL rendering is required for transparency (at least on GNU/Linux)

plot3c = @(varargin)plot3(varargin{1},real(varargin{2}),imag(varargin{2}),varargin{3:end});

ax(1) = subplot(221);
% plot the signal
plot3c(absc,x(tt))
grid
axis(A)
set(gca,Axprops,Axpropvals);
subplot_label('a')

ax(2) = subplot(222);
% plot the signal
plot3c(absc,x(tt))
% and the envelope curves
for k = 1:2*Ndirs1
    plot3c(absc,env1(k,tt),'k','LineWidth',L)
end
set(gca,Axprops,Axpropvals);
grid
% plot the tube envelope
h = surf(PX,real(PYZ),imag(PYZ),gray,'VertexNormals',VN);
set(h,propNames,propVals);
axis(A)
l1 = light('Position',[25,15,20]);
l2 = light('Position',[-50,10,-18]);
subplot_label('b')

ax(3) = subplot(223);
% plot the mean of the envelope curves
plot3c(absc,x(tt)-moy(tt))
% and the zero axis
plot3c(absc,zeros(1,length(absc)),'k')
axis(A)
set(gca,Axprops,Axpropvals);
xlabel('Time')
grid
subplot_label('c')

ax(4) = subplot(224);
plot3c(absc,moy(tt))
axis(A)
set(gca,Axprops,Axpropvals);
xlabel('Time')
grid
subplot_label('d')

hlink = linkprop(ax,{'CameraPosition','CameraUpVector'});
