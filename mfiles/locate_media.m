function mediainfo = locate_media(parfnm,varargin)

% locate media index in mpi threads
% Author:   Yuanhang Huo
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.06.06

gtstart  = 1;
gtcount  = -1;
gtstride = 1;

%-- flags --
n=1;
while n<=nargin-1
    
    if numel(varargin{n})==1 | ~isnumeric(varargin{n})
        switch varargin{n}
            case 'start'
                gsubs=varargin{n+1}; n=n+1;
                if length(gsubs)==4
                    gtstart=gsubs(4);gsubs=gsubs(1:3);
                end
            case 'count'
                gsubc=varargin{n+1}; n=n+1;
                if length(gsubc)==4
                    gtcount=gsubc(4);gsubc=gsubc(1:3);
                end
            case 'stride'
                gsubt=varargin{n+1}; n=n+1;
                if length(gsubt)==4
                    gtstride=gsubt(4);gsubt=gsubt(1:3);
                end
            case 'mediadir'
                media_dir=varargin{n+1}; n=n+1;
        end
    end
    
    n=n+1;
    
end

% check parameter file exist
if ~ exist(parfnm,'file')
    error([mfilename ': file ' parfnm ' does not exist']);
end

% read parameters file
par=loadjson(parfnm);
media_subs=[1, 1, 1];
media_subc=[-1,-1,-1];
media_subt=[1, 1, 1];
snap_tinv=1;
ngijk=[par.number_of_total_grid_points_x,...
       par.number_of_total_grid_points_y,...
       par.number_of_total_grid_points_z];

% reset count=-1 to total number
indx=find(media_subc==-1);
media_subc(indx)=fix((ngijk(indx)-media_subs(indx))./media_subt(indx))+1;
snap_sube=media_subs+(media_subc-1).*media_subt;

indx=find(gsubc==-1);
gsubc(indx)=fix((media_subc(indx)-gsubs(indx))./gsubt(indx))+1;
gsube=gsubs+(gsubc-1).*gsubt;

% search the nc file headers to locate the threads/processors
mediaprefix='media';
medialist=dir([media_dir,'/',mediaprefix,'*.nc']);
n=1;
for i=1:length(medialist)
    
    medianm=[media_dir,'/',medialist(i).name];
    xyzs=double(nc_attget(medianm,nc_global,'global_index_of_first_physical_points'));
    xs=xyzs(1);
    ys=xyzs(2);
    xyzc=double(nc_attget(medianm,nc_global,'count_of_physical_points'));
    xc=xyzc(1);
    yc=xyzc(2);
    xarray=[xs:xs+xc-1];
    yarray=[ys:ys+yc-1];
    if length(find(xarray>=gsubs(1)-1 & xarray<=gsube(1)-1)) ~= 0 && ...
       length(find(yarray>=gsubs(2)-1 & yarray<=gsube(2)-1)) ~= 0
        px(n)=str2num(medialist(i).name( strfind(medialist(i).name,'px' )+2 : ...
                                         strfind(medialist(i).name,'_py')-1));
        py(n)=str2num(medialist(i).name( strfind(medialist(i).name,'py' )+2 : ...
                                         strfind(medialist(i).name,'.nc')-1));
        n=n+1;
    end
    
end

% retrieve the snapshot information
nthd=0;
for ip=1:length(px)
    
    nthd=nthd+1;
    
    medianm=[media_dir,'/',mediaprefix,'_px',num2str(px(ip)),...
            '_py',num2str(py(ip)),'.nc'];
    xyzs=double(nc_attget(medianm,nc_global,'global_index_of_first_physical_points'));
    xs=xyzs(1);
    ys=xyzs(2);
    zs=xyzs(3);
    xyzc=double(nc_attget(medianm,nc_global,'count_of_physical_points'));
    xc=xyzc(1);
    yc=xyzc(2);
    zc=xyzc(3);
    xe=xs+xc-1;
    ye=ys+yc-1;
    ze=zs+zc-1;
    
    gxarray=gsubs(1):gsubt(1):gsube(1);
    gxarray=gxarray-1;
    gyarray=gsubs(2):gsubt(2):gsube(2);
    gyarray=gyarray-1;
    gzarray=gsubs(3):gsubt(3):gsube(3);
    gzarray=gzarray-1;
    
    mediainfo(nthd).thisid=[px(ip),py(ip)];
    i=find(gxarray>=xs & gxarray<=xe);
    j=find(gyarray>=ys & gyarray<=ye);
    k=find(gzarray>=zs & gzarray<=ze);
    mediainfo(nthd).indxs=[i(1),j(1),k(1)];
    mediainfo(nthd).indxe=[i(end),j(end),k(end)];
    mediainfo(nthd).indxc=mediainfo(nthd).indxe-mediainfo(nthd).indxs+1;
    
    mediainfo(nthd).subs=[ gxarray(i(1))-xs+1, ...
                           gyarray(j(1))-ys+1, ...
                           gzarray(k(1))-zs+1 ];
    mediainfo(nthd).subc=mediainfo(nthd).indxc;
    mediainfo(nthd).subt=gsubt;
    
    mediainfo(nthd).wsubs=double(nc_attget(medianm,nc_global,'local_index_of_first_physical_points'))...
        +(mediainfo(nthd).subs-1).*media_subt+1;
    mediainfo(nthd).wsubc=mediainfo(nthd).indxc;
    mediainfo(nthd).wsubt=media_subt.*gsubt;
    
    mediainfo(nthd).tinv=snap_tinv;
    
    mediainfo(nthd).fnmprefix=mediaprefix;
    
    mediainfo(nthd).ttriple=[gtstart,gtcount,gtstride];
    
end


end

