function [x,y,z] = gather_coord(snapinfo,varargin)

% gather coordinates of the snapshot
% Author:   Yuanhang Huo
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.05.31

%-- flags --
nargs=nargin-1;
n=1;
while n<=nargs
    
    if numel(varargin{n})==1 | ~isnumeric(varargin{n})
        switch varargin{n}
            case 'coorddir'
                output_dir=varargin{n+1}; n=n+1;
        end
    end
    
    n=n+1;
    
end

% check
if ~ exist(output_dir,'dir')
    error([mfilename ': directory ' output_dir ' does not exist']);
end

% load
coordprefix='coord';
nthd=length(snapinfo);
for n=1:nthd
    
    n_i=snapinfo(n).thisid(1); n_j=snapinfo(n).thisid(2);
    i1=snapinfo(n).indxs(1); j1=snapinfo(n).indxs(2); k1=snapinfo(n).indxs(3);
    i2=snapinfo(n).indxe(1); j2=snapinfo(n).indxe(2); k2=snapinfo(n).indxe(3);
    subs=snapinfo(n).wsubs;
    subc=snapinfo(n).wsubc;
    subt=snapinfo(n).wsubt;
    fnm_coord=[output_dir,'/',coordprefix,'_px',num2str(n_i),'_py',num2str(n_j),'.nc'];
    
    if ~ exist(fnm_coord,'file')
       error([mfilename ': file ' fnm_coord 'does not exist']);
    end

    subs=fliplr(subs);subc=fliplr(subc);subt=fliplr(subt);
    x(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'x',subs-1,subc,subt);
    y(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'y',subs-1,subc,subt);
    z(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'z',subs-1,subc,subt);
    
end

x=double(permute(x,[3 2 1]));
y=double(permute(y,[3 2 1]));
z=double(permute(z,[3 2 1]));


end