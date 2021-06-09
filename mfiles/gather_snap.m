function [V,varargout] =gather_snap(snapinfo,nlayer,varnm,varargin)

% gather snapshot data
% Author:   Yuanhang Huo
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.05.31

%-- flags --
nargs=nargin-3;
n=1;

while n<=nargs
    if numel(varargin{n})==1 | ~isnumeric(varargin{n})
        switch varargin{n}
            case 'snapdir'
                snap_dir=varargin{n+1}; n=n+1;
        end
    end   
    n=n+1;
end

% check dir exists
if ~ exist(snap_dir,'dir')
    error([mfilename ': directory ' snap_dir ' does not exist']);
end

% load
nthd=length(snapinfo);
for n=1:nthd
    
    n_i=snapinfo(n).thisid(1); n_j=snapinfo(n).thisid(2);
    i1=snapinfo(n).indxs(1); j1=snapinfo(n).indxs(2); k1=snapinfo(n).indxs(3);
    i2=snapinfo(n).indxe(1); j2=snapinfo(n).indxe(2); k2=snapinfo(n).indxe(3);
    subs=snapinfo(n).subs;
    subc=snapinfo(n).subc;
    subt=snapinfo(n).subt;
    fnm_snap=[snap_dir,'/',snapinfo(n).fnmprefix,'_px',num2str(n_i),...
              '_py',num2str(n_j),'.nc'];
    if ~ exist(fnm_snap)
       error([mfilename ': file ',fnm_snap, ' does not exist']);
    end
    tdim=nc_getdiminfo(fnm_snap,'time');
    if tdim.Length==0 | (nlayer-1)-1>=tdim.Length
       error([num2str(nlayer) 'th layer is beyond current time dim (' ...
            num2str(tdim.Length) ') in ' fnm_snap]);
    end
    % get data
    V(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_snap,varnm, ...
          [nlayer-1,fliplr(subs)-1],[1,fliplr(subc)],[1,fliplr(subt)]);
    t=nc_varget(fnm_snap,'time',[nlayer-1],[1]);

end

V=double(permute(V,[3,2,1]));

if nargout>=2
    varargout(1)={t};
end
if nargout>=3
    varargout(2)={varnm};
end


end
