function sid=vef_surf_light(h,x,y,z,v,cMap,varargin)

% ----------------------- parameter -----------------------
args=varargin; nargs=numel(args);
n=1;
while n<=nargs

if ~ischar(args{n}), n=n+1; continue; end

switch args{n}
case 'bgcolor', bgcolor=args{n+1};
case 'daspect', scl_daspect=args{n+1};
end

n=n+1;
end

% ----------------------- surf -----------------------
sid=surf(x,y,z,v);
set(sid,'DiffuseStrength',1.0,'SpecularStrength',0.2, ...
    'SpecularExponent',50, ...
    'SpecularColorReflectance',0.1)

% ----------------------- annotation -----------------------
if exist('bgcolor','var')
   set(h,'color',bgcolor);
   set(gca,'color',bgcolor);
end

if exist('scl_daspect'); daspect(scl_daspect); end

axis tight

% -- colormap --
colormap(cMap);
%vef_colormap(@jet);

% -- shading --
%shading interp;
shading flat;

%vef_light;

% -- light --
view(-20,35);
set(gca,'box','off');
%if ~exist('lid','var'),
   lid=camlight('headlight','local');
%end
camlight(lid,70,8,'local')
%camlight(lid,-30,60)
lighting phong
set(gca,'box','on');
%lid=camlight('headlight','infinite')
%camlight(lid,40,15,'local')
%camlight(lid,30,40) 

% -- view angle --
view(0,90)
