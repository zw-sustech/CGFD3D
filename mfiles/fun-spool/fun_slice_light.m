function [h,sid]=fun_slice_light(xd,yd,zd,varargin)

global X Y Z V

% ----------------------- parameter -----------------------
flag_jetwr = 0;

args=varargin; nargs=numel(args);
n=1;
while n<=nargs

if ~ischar(args{n}), n=n+1; continue; end

switch args{n}
case 'bgcolor', bgcolor=args{n+1};
case 'caxis', scl_caxis=args{n+1};
case 'xlim', scl_xlim=args{n+1};
case 'ylim', scl_ylim=args{n+1};
case 'zlim', scl_zlim=args{n+1};
case 'daspect', scl_daspect=args{n+1};
case 'jetwr', flag_jetwr=1; 
case 'station' 
  STLO=args{n+1}(1);STLA=args{n+1}(2);
case 'event' 
  EVLO=args{n+1}(1);EVLA=args{n+1}(2);EVDP=args{n+1}(3);
end

n=n+1;
end

% ----------------------- slice -----------------------
h=figure;
set(h,'renderer','zbuffer');
set(h,'menubar','none');
set(h,'toolbar','figure');
%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf,'PaperUnits','points');
%set(gcf,'PaperPosition',[0 0 1024 768]);
%set(0, 'DefaultFigurePaperType', 'A4');

nslc=numel(xd);

for m=1:nslc
    sid(m)=slice(X,Y,Z,V,xd{m},yd{m},zd{m});
    %sid(m)=slice(Y,X,Z,V,xd{m},yd{m},zd{m});
    hold on
end

% ----------------------- annotation -----------------------
%set(sid(1:end-1),'DiffuseStrength',1.0,'SpecularStrength',0.2, ...
set(sid(1:end),'DiffuseStrength',1.0,'SpecularStrength',0.2, ...
    'SpecularExponent',50, ...
    'SpecularColorReflectance',0.1)
%set(sid(3),'DiffuseStrength',1.0,'SpecularStrength',0.2, ...
%    'SpecularExponent',50, ...
%    'SpecularColorReflectance',0.1)

if exist('bgcolor','var')
   set(h,'color',bgcolor);
   set(gca,'color',bgcolor);
end

if exist('scl_daspect'); daspect(scl_daspect); end

grid on
box off
axis tight

%axis image
%shading flat;
shading interp;
if flag_jetwr==1
   colormap('jetwr');
end
if exist('scl_caxis','var'), caxis(scl_caxis); end
if exist('scl_xlim','var'), xlim(scl_xlim); end
if exist('scl_ylim','var'), ylim(scl_ylim); end
if exist('scl_zlim','var'), zlim(scl_zlim); end
set(gca,'zdir','reverse');

%pid(1)=plot3(EVLO,EVLA,20,'o','Markersize',8,'MarkerEdgeColor','k','MarkerFaceColor','w')
if exist('EVLO','var')
sid(end+1)=plot3(EVLO,EVLA,EVDP,'p','Markersize',8,'MarkerEdgeColor','k','MarkerFaceColor','w');
end
if exist('STLO','var')
sid(end+1)=plot3(STLO,STLA,0,'v','Markersize',8,'MarkerEdgeColor','k','MarkerFaceColor','w');
end


% -- light effect ---
%view(-20,30);
view(10,30);
%view(20,30);
set(gca,'box','off');
%if ~exist('lid','var'),
   %lid=camlight('headlight','local');
   lid=camlight('headlight');
%end
%camlight(lid,70,8,'local')
%camlight(lid,50,-8)
%c1amlight(lid,10,20)
camlight(lid,0,10);
lighting phong
%camlight(lid,15,15)
%set(gca,'box','on');
%lid=camlight('headlight','infinite')
%camlight(lid,40,15,'local')
%camlight(lid,30,40) 

view(20,30);

