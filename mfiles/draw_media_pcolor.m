% Draw a cross-section of media by using pcolor style
% Author:       Yuanhang Huo
% Email:        yhhuo@mail.ustc.edu.cn
% Affiliation:  University of Science and Technology of China
% Date:         2021.05.31

% file and path name
parfnm='./project/test.json';
output_dir='./project/output';

% which snapshot to plot
id=1;
subs=[1,1,50];
subc=[-1,-1,1];
subt=[1,1,1];

% variable and time to plot
varnm='Vp';

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_clb    = 1;
scl_daspect =[1 1 1];
clrmp       = 'parula';

% load snapshot data
snapinfo=locate_media(parfnm,'start',subs,'count',subc,'stride',subt,'outdir',output_dir);
% get coordinate data
[x,y,z]=gather_coord(snapinfo,'coorddir',output_dir);
nx=size(x,1);
ny=size(x,2);
nz=size(x,3);
% coordinate unit
str_unit='m';
if flag_km
   x=x/1e3;
   y=y/1e3;
   z=z/1e3;
   str_unit='km';
end

% load media data
switch varnm
case 'Vp'
   rho=gather_media(snapinfo,'rho','mediadir',dir_media);
   mu=gather_media(snapinfo,'mu','mediadir',dir_media);
   lambda=gather_media(snapinfo,'lambda','mediadir',dir_media);
   v=( (lambda+2*mu)./rho ).^0.5;
   v=v/1e3;
case 'Vs'
   rho=gather_media(snapinfo,'rho','mediadir',dir_media);
   mu=gather_media(snapinfo,'mu','mediadir',dir_media);
   v=( mu./rho ).^0.5;
   v=v/1e3;
case 'rho'
   v=gather_media(snapinfo,varnm,'mediadir',dir_media);
   v=v/1e3;
otherwise
   v=gather_media(snapinfo,varnm,'mediadir',dir_media);
end


% figure plot
hid=figure;
set(hid,'BackingStore','on');
set(hid,'renderer','painters');
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperUnits','points');
set(gcf,'PaperPosition',[0 0 800 800]);

if nx==1
   if flag_emlast
      sid=pcolor((flipud(permute(squeeze(y),[2 1]))), ...
                 (flipud(permute(squeeze(z),[2 1]))), ...
                 (flipud(permute(squeeze(v),[2 1]))));
   else
      sid=pcolor(permute(squeeze(y),[2 1]), ...
                 permute(squeeze(z),[2 1]), ...
                 permute(squeeze(v),[2 1]));
   end
   xlabel(['Y axis (' str_unit ')']);
   ylabel(['Z axis (' str_unit ')']);
   
elseif ny==1
   if flag_emlast
      %v(:,:,1)=0; v(:,:,end)=0;
      sid=pcolor(flipud(permute(squeeze(x),[2 1])), ...
                 flipud(permute(squeeze(z),[2 1])), ...
                 flipud(permute(squeeze(v),[2 1])));
   else
      sid=pcolor(permute(squeeze(x),[2 1]), ...
                 permute(squeeze(z),[2 1]), ...
                 permute(squeeze(v),[2 1]));
   end
   xlabel(['X axis (' str_unit ')']);
   ylabel(['Z axis (' str_unit ')']);
   %set(gca,'ydir','reverse');
   
else
   if flag_emlast
      sid=pcolor(flipud(permute(squeeze(x),[2 1])), ...
                 flipud(permute(squeeze(y),[2 1])), ...
                 flipud(permute(squeeze(v),[2 1])));
   else
      sid=pcolor(permute(squeeze(x),[2 1]), ...
                 permute(squeeze(y),[2 1]), ...
                 permute(squeeze(v),[2 1]));
   end
   xlabel(['X axis (' str_unit ')']);
   ylabel(['Y axis (' str_unit ')']);
end

% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis tight

% colormap and colorbar
if exist('clrmp')
    colormap(jetwr);
end
if exist('scl_caxis','var')
    caxis(scl_caxis);
end
if flag_clb
    cid=colorbar;
end

% shading
%shading interp;
shading flat;

% title
if flag_title
    title(varnm);
end

% save figure
if flag_print==1
   print(gcf,[varnm '.png'],'-dpng');
end











