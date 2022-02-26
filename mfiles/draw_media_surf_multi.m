% Draw multi cross-sections of media by using surf style
% Author:       Yuanhang Huo
% Email:        yhhuo@mail.ustc.edu.cn
% Affiliation:  University of Science and Technology of China
% Date:         2021.06.09

clear all;
addmypath
% -------------------------- parameters input -------------------------- %
% file and path name
%media_type = 'ac_iso';
media_type = 'el_iso';
parfnm='../project/test.json'
output_dir='../project/output'

% media profiles to plot
% profile 1
subs{1}=[50,1,1];      % start from index '1'
subc{1}=[1,-1,-1];     % '-1' to plot all points in this dimension
subt{1}=[1,1,1];
% profile 2
subs{2}=[1,1,50];      % start from index '1'
subc{2}=[-1,-1,1];     % '-1' to plot all points in this dimension
subt{2}=[1,1,1];
% profile 3
subs{3}=[1,50,1];      % start from index '1'
subc{3}=[-1,1,-1];     % '-1' to plot all points in this dimension
subt{3}=[1,1,1];

% variable to plot
% 'Vp', 'Vs', 'rho', 'lambda', 'mu'
varnm='Vp';

% figure control parameters
flag_km     = 1;
flag_emlast = 1;
flag_print  = 0;
flag_clb    = 1;
flag_title  = 1;
scl_daspect =[1 1 1];
clrmp       = 'parula';
% ---------------------------------------------------------------------- %



% figure plot
hid=figure;
set(hid,'BackingStore','on');

% load data and plot
for i=1:length(subs)
    
    % locate media
    mediainfo{i}=locate_media(parfnm,'start',subs{i},'count',subc{i},'stride',subt{i},'mediadir',output_dir);
    
    % get coordinate data
    [x{i},y{i},z{i}]=gather_coord(mediainfo{i},'coorddir',output_dir);
    nx{i}=size(x{i},1);
    ny{i}=size(x{i},2);
    nz{i}=size(x{i},3);
    
    % coordinate unit
    str_unit='m';
    if flag_km
        x{i}=x{i}/1e3;
        y{i}=y{i}/1e3;
        z{i}=z{i}/1e3;
        str_unit='km';
    end
    
    % gather media
    switch varnm
        case 'Vp'
            rho=gather_media(mediainfo{i},'rho',output_dir);
            if strcmp(media_type,'ac_iso') == 1
               kappa=gather_media(mediainfo{i},'kappa',output_dir);
               v{i}=( kappa ./rho ).^0.5;
            elseif strcmp(media_type,'el_iso') == 1
               mu=gather_media(mediainfo{i},'mu',output_dir);
               lambda=gather_media(mediainfo{i},'lambda',output_dir);
               v{i}=( (lambda+2*mu)./rho ).^0.5;
            end
            v{i}=v{i}/1e3;
        case 'Vs'
            rho=gather_media(mediainfo{i},'rho',output_dir);
            mu=gather_media(mediainfo{i},'mu',output_dir);
            v{i}=( mu./rho ).^0.5;
            v{i}=v{i}/1e3;
        case 'rho'
            v{i}=gather_media(mediainfo{i},varnm,output_dir);
            v{i}=v{i}/1e3;
        otherwise
            v{i}=gather_media(mediainfo{i},varnm,output_dir);
    end
    
    % media show
    if flag_emlast
        sid{i}=surf(squeeze(permute(x{i},[2 1 3])), ...
            squeeze(permute(y{i},[2 1 3])), ...
            squeeze(permute(z{i},[2 1 3])), ...
            squeeze(permute(v{i},[2 1 3])));
    else
        sid{i}=surf(flipdim(squeeze(permute(x{i},[2 1 3])),3), ...
            flipdim(squeeze(permute(y{i},[2 1 3])),3), ...
            flipdim(squeeze(permute(z{i},[2 1 3])),3), ...
            flipdim(squeeze(permute(v{i},[2 1 3])),3));
    end

    hold on;
    
end

xlabel(['X axis (' str_unit ')']);
ylabel(['Y axis (' str_unit ')']);
zlabel(['Z axis (' str_unit ')']);

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


