%bivariate_EMD_mean_definitions.m 
%shows how the center of the tube envelope is defined for each algorithm
%reproduces Fig. 2 in "Bivariate Empirical Mode Decomposition", G. Rilling,
%P. Flandrin, P. Goncalves and J. M. Lilly, IEEE Signal Processing Letters
%
%Note that the script allows to easily change the signal model, number of
%directions and what is displayed
%
%G. Rilling 3/2007 email:  gabriel.rilling@ens-lyon.fr

%number of directions
Ndirs = 4;

%display options
show_errorbars = 1;
show_tangents = 1;
show_midtangents = 1;
show_polygon = 1;

i = sqrt(-1); % just in case...

% example 1: elliptical signal
signal = @(t)dirstretch(exp(i*t),pi/3,2);

% example 2: less symetrical signal
% signal = @(t)exp(i*pi/5)*(exp(i*t)+.1*exp(-2*i*t)+.1*exp(-i*t));

% NOTE: you can try your own examples provided the cycle shape is convex


x = signal(absc);
absc = linspace(0,2*pi,1e3);
dirs = (0:Ndirs-1)*2*pi/Ndirs;
% get indices of extreme points
inds = [];
for cdi=1:length(dirs) 
  dir = dirs(cdi);
  [val,ind] = max(real(exp(-i*dir)*x));
  inds(cdi) = ind(1);
end

%display parameters
L = 2;
M = 5;
Length = 10;
A = 2.5*[-1,1,-1,1];
label_pos = [.05,.05];
subplot_label = @(label)text(A(1)+(A(2)-A(1))*label_pos(1),A(3)+(A(4)-A(3))*(1-label_pos(2)),['(',label,')']);
small_bar = .25;
big_bar = .7;
mid_bar = (big_bar+small_bar)/2;

% draws schematic errorbar
plotbar = @(x,dirs,sb,bb) plot([x-bb*dirs,x-sb*dirs*i;x+bb*dirs,x+sb*dirs*i],'Color','k','LineWidth',L);
if ~show_errorbars
  plotbar = @(x,dirs,sb,bb)[];
end
% draws bigpoints
plotpoint = @(x,col) plot(x,[col,'o'],'MarkerSize',M,'MarkerFaceColor',col);
% draws a line in direction dir passing through point x 
tline = @(x,dir,varargin) plot(real([x-Length*exp(i*dir);x+Length*exp(i*dir)]),imag([x-Length*exp(i*dir);x+Length*exp(i*dir)]),varargin{:});

figtitle = 'Definition of the mean of the envelope';
figure('name',figtitle)

% first algorithm
subplot(121)
set(gca,'NextPlot','add','Box','on','XTick',[],'YTick',[])
plot(x,'k')
plotbar(x(inds),exp(i*dirs),big_bar,small_bar);
plotpoint(x(inds),'b')
if show_polygon
  plot(x([inds,inds(1)]),'b')
end
plotbar(mean(x(inds)),1,mid_bar,mid_bar);
plotpoint(mean(x(inds)),'r')
subplot_label('a')
axis(A);
axis square
title('algorithm 1')

%second algorithm
subplot(122)
set(gca,'NextPlot','add','Box','on','XTick',[],'YTick',[])
plot(x,'k')
if show_tangents
  tline(x(inds),dirs+pi/2,'b')
end
c = exp(i*dirs).*(real(exp(-i*dirs).*x(inds))); % <- quantity used to compute the center of the tube
plotbar(x(inds),exp(i*dirs),0,small_bar);
plotpoint(x(inds),'b')
% draw lines halfway between tangents if Ndirs is even
if mod(Ndirs,2)==0 && show_midtangents
  inds2 = 1:Ndirs/2;
  tline(exp(i*dirs(inds2)).*(real(exp(-i*dirs(inds2)).*(x(inds(inds2))+x(inds(inds2+Ndirs/2)))/2)),dirs(inds2)+pi/2,'--b')
end
% draw the center
center = 2*mean(c);
plotbar(center,1,small_bar,small_bar);
plotpoint(center,'r')
subplot_label('b')
axis(A);
axis square
title('algorithm 2')