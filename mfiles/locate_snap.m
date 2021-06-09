function snapinfo = locate_snap(parfnm,id,varargin)

% locate snapshot index in mpi threads
% Author:   Yuanhang Huo
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.05.31

gtstart  = 1;
gtcount  = -1;
gtstride = 1;

%-- flags --
n=1;
while n<=nargin-2
    
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
            case 'snapdir'
                snap_dir=varargin{n+1}; n=n+1;
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
snap_subs=par.snapshot{id}.grid_index_start+1;
snap_subc=par.snapshot{id}.grid_index_count;
snap_subt=par.snapshot{id}.grid_index_incre;
snap_tinv=par.snapshot{id}.time_index_incre;
ngijk=[par.number_of_total_grid_points_x,...
       par.number_of_total_grid_points_y,...
       par.number_of_total_grid_points_z];

% reset count=-1 to total number
indx=find(snap_subc==-1);
snap_subc(indx)=fix((ngijk(indx)-snap_subs(indx))./snap_subt(indx))+1;
snap_sube=snap_subs+(snap_subc-1).*snap_subt;

indx=find(gsubc==-1);
gsubc(indx)=fix((snap_subc(indx)-gsubs(indx))./gsubt(indx))+1;
gsube=gsubs+(gsubc-1).*gsubt;

% search the nc file headers to locate the threads/processors
snapprefix=par.snapshot{id}.name;
snaplist=dir([snap_dir,'/',snapprefix,'*.nc']);
n=1;
for i=1:length(snaplist)
    
    snapnm=[snap_dir,'/',snaplist(i).name];
    xyzs=double(nc_attget(snapnm,nc_global,'first_index_to_snapshot_output'));
    xs=xyzs(1);
    ys=xyzs(2);
    xc=nc_getdiminfo(snapnm,'i');
    xc=xc.Length;
    yc=nc_getdiminfo(snapnm,'j');
    yc=yc.Length;
    xarray=[xs:xs+xc-1];
    yarray=[ys:ys+yc-1];
    if length(find(xarray>=gsubs(1)-1 & xarray<=gsube(1)-1)) ~= 0 && ...
       length(find(yarray>=gsubs(2)-1 & yarray<=gsube(2)-1)) ~= 0
        px(n)=str2num(snaplist(i).name( strfind(snaplist(i).name,'px' )+2 : ...
                                        strfind(snaplist(i).name,'_py')-1));
        py(n)=str2num(snaplist(i).name( strfind(snaplist(i).name,'py' )+2 : ...
                                        strfind(snaplist(i).name,'.nc')-1));
        n=n+1;
    end
    
end

% retrieve the snapshot information
nthd=0;
for ip=1:length(px)
    
    nthd=nthd+1;
    
    snapnm=[snap_dir,'/',snapprefix,'_px',num2str(px(ip)),...
            '_py',num2str(py(ip)),'.nc'];
    xyzs=double(nc_attget(snapnm,nc_global,'first_index_to_snapshot_output'));
    xs=xyzs(1);
    ys=xyzs(2);
    zs=xyzs(3);
    xc=nc_getdiminfo(snapnm,'i');
    xc=xc.Length;
    yc=nc_getdiminfo(snapnm,'j');
    yc=yc.Length;
    zc=nc_getdiminfo(snapnm,'k');
    zc=zc.Length;
    xe=xs+xc-1;
    ye=ys+yc-1;
    ze=zs+zc-1;
    
    gxarray=gsubs(1):gsubt(1):gsube(1);
    gxarray=gxarray-1;
    gyarray=gsubs(2):gsubt(2):gsube(2);
    gyarray=gyarray-1;
    gzarray=gsubs(3):gsubt(3):gsube(3);
    gzarray=gzarray-1;
    
    snapinfo(nthd).thisid=[px(ip),py(ip)];
    i=find(gxarray>=xs & gxarray<=xe);
    j=find(gyarray>=ys & gyarray<=ye);
    k=find(gzarray>=zs & gzarray<=ze);
    snapinfo(nthd).indxs=[i(1),j(1),k(1)];
    snapinfo(nthd).indxe=[i(end),j(end),k(end)];
    snapinfo(nthd).indxc=snapinfo(nthd).indxe-snapinfo(nthd).indxs+1;
    
    snapinfo(nthd).subs=[ gxarray(i(1))-xs+1, ...
                          gyarray(j(1))-ys+1, ...
                          gzarray(k(1))-zs+1 ];
    snapinfo(nthd).subc=snapinfo(nthd).indxc;
    snapinfo(nthd).subt=gsubt;
    
    snapinfo(nthd).wsubs=double(nc_attget(snapnm,nc_global,'first_index_in_this_thread_with_ghosts'))...
        +(snapinfo(nthd).subs-1).*snap_subt+1;
    snapinfo(nthd).wsubc=snapinfo(nthd).indxc;
    snapinfo(nthd).wsubt=snap_subt.*gsubt;
    
    snapinfo(nthd).tinv=snap_tinv;
    
    snapinfo(nthd).fnmprefix=snapprefix;
    
    snapinfo(nthd).ttriple=[gtstart,gtcount,gtstride];
    
end


end

