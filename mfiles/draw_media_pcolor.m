% Draw a cross-section of media by using pcolor style
% Author:       Yuanhang Huo
% Email:        yhhuo@mail.ustc.edu.cn
% Affiliation:  University of Science and Technology of China
% Date:         2021.06.06

clear all;

% -------------------------- parameters input -------------------------- %
% file and path name
%media_type = 'ac_iso';
%parfnm='/home/zhangw/work/cgfd3d-wave-ac/02/test.json'
%output_dir='/home/zhangw/work/cgfd3d-wave-ac/02/output'
media_type = 'el_iso';
parfnm='/home/zhangw/work/cgfd3d-wave-el/04/test.json'
output_dir='/home/zhangw/work/cgfd3d-wave-el/04/output'

%media_type = 'el_vti';

% which media profile to plot
subs=[1,50,1];      % start from index '1'
subc=[-1,1,-1];     % '-1' to plot all points in this dimension
subt=[1,1,1];

%subs=[50,1,1];      % start from index '1'
%subc=[1,-1,-1];     % '-1' to plot all points in this dimension
%subt=[1,1,1];

% variable to plot
% 'Vp', 'Vs', 'rho', 'lambda', 'mu'
varnm='Vp';

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_clb    = 1;
flag_title  = 1;
scl_daspect = [1 1 1];
clrmp       = 'parula';
% ---------------------------------------------------------------------- %



% load media data
mediainfo=locate_media(parfnm,'start',subs,'count',subc,'stride',subt,'mediadir',output_dir);
% get coordinate data
[x,y,z]=gather_coord(mediainfo,'coorddir',output_dir);
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
        rho=gather_media(mediainfo,'rho',output_dir);
        if strcmp(media_type,'ac_iso') == 1
          kappa=gather_media(mediainfo,'kappa',output_dir);
          v=( (kappa)./rho ).^0.5;
        elseif strcmp(media_type,'el_iso') == 1
          mu=gather_media(mediainfo,'mu',output_dir);
          lambda=gather_media(mediainfo,'lambda',output_dir);
          v=( (lambda+2*mu)./rho ).^0.5;
        end
        v=v/1e3;
    case 'Vs'
        rho=gather_media(mediainfo,'rho',output_dir);
        mu=gather_media(mediainfo,'mu',output_dir);
        v=( mu./rho ).^0.5;
        v=v/1e3;
    case 'rho'
        v=gather_media(mediainfo,varnm,output_dir);
        v=v/1e3;
    otherwise
        v=gather_media(mediainfo,varnm,output_dir);
end


% figure plot
hid=figure;
set(hid,'BackingStore','on');

% media show
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

set(gca,'layer','top');
set(gcf,'color','white','renderer','painters');

% shading
% shading interp;
shading flat;
% colorbar range/scale
if exist('scl_caxis','var')
    caxis(scl_caxis);
end
% axis daspect
if exist('scl_daspect')
    daspect(scl_daspect);
end
axis tight
% colormap and colorbar
if exist('clrmp')
    colormap(clrmp);
end
if flag_clb
    cid=colorbar;
    if strcmp(varnm,'Vp') || strcmp(varnm,'Vs')
        cid.Label.String='(km/s)';
    end
    if strcmp(varnm,'rho')
        cid.Label.String='g/cm^3';
    end
end

% title
if flag_title
    title(varnm);
end

% save and print figure
if flag_print
    width= 500;
    height=500;
    set(gcf,'paperpositionmode','manual');
    set(gcf,'paperunits','points');
    set(gcf,'papersize',[width,height]);
    set(gcf,'paperposition',[0,0,width,height]);
    print(gcf,[varnm '.png'],'-dpng');
end


