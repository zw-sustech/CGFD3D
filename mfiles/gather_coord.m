function [x,y,z] = gather_coord(coordinfo,varargin)

% gather coordinates of the snapshot
% Author:   Yuanhang Huo
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.05.31

% ChangeLog:
% 20210814: read both cart and curv coord, --zw

%-- flags --
nargs=nargin-1;
n=1;
while n<=nargs
    
    if numel(varargin{n})==1 | ~isnumeric(varargin{n})
        switch varargin{n}
            case 'coorddir'
                coord_dir=varargin{n+1}; n=n+1;
        end
    end
    
    n=n+1;
    
end

% check
if ~ exist(coord_dir,'dir')
    error([mfilename ': directory ' coord_dir ' does not exist']);
end

% load
coordprefix='coord';
nthd=length(coordinfo);
for n=1:nthd
    
    n_i=coordinfo(n).thisid(1); n_j=coordinfo(n).thisid(2);
    i1=coordinfo(n).indxs(1); j1=coordinfo(n).indxs(2); k1=coordinfo(n).indxs(3);
    i2=coordinfo(n).indxe(1); j2=coordinfo(n).indxe(2); k2=coordinfo(n).indxe(3);
    subs=coordinfo(n).wsubs;
    subc=coordinfo(n).wsubc;
    subt=coordinfo(n).wsubt;
    fnm_coord=[coord_dir,'/',coordprefix,'_px',num2str(n_i),'_py',num2str(n_j),'.nc'];
    
    if ~ exist(fnm_coord,'file')
       error([mfilename ': file ' fnm_coord 'does not exist']);
    end

    subs=fliplr(subs);subc=fliplr(subc);subt=fliplr(subt);

    % check dimension size of x,y,z to detmine cart or curv
    xvar_info = ncinfo(fnm_coord,'x');
    num_dim_x = length(xvar_info.Dimensions);
    if num_dim_x == 1 % cart grid
      x1d = nc_varget(fnm_coord,'x',subs(3)-1,subc(3),subt(3));
      x3d = repmat(x1d,1,(j2-j1+1),(k2-k1+1));
      x(k1:k2,j1:j2,i1:i2) = permute(x3d,[3 2 1]);
    else % curv grid
      x(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'x',subs-1,subc,subt);
    end

    %- y coord
    yvar_info = ncinfo(fnm_coord,'y');
    num_dim_y = length(yvar_info.Dimensions);
    if num_dim_y == 1 % cart grid
      y1d = nc_varget(fnm_coord,'y',subs(2)-1,subc(2),subt(2));
      y3d = repmat(y1d,1,(i2-i1+1),(k2-k1+1)); % y1d is a col, can't be rep as 1st dim
      y(k1:k2,j1:j2,i1:i2) = permute(y3d,[3 1 2]);
    else % curv grid
      y(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'y',subs-1,subc,subt);
    end

    %- z coord
    zvar_info = ncinfo(fnm_coord,'z');
    num_dim_z = length(zvar_info.Dimensions);
    if num_dim_z == 1 % cart grid
      z1d = nc_varget(fnm_coord,'z',subs(1)-1,subc(1),subt(1));
      z3d = repmat(z1d,1,(j2-j1+1),(i2-i1+1));
      %z(k1:k2,j1:j2,i1:i2) = permute(z3d,[1 2 3]);
      z(k1:k2,j1:j2,i1:i2) = z3d;
    else % curv grid
      z(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'z',subs-1,subc,subt);
    end

    % only for curv
    %x(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'x',subs-1,subc,subt);
    %y(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'y',subs-1,subc,subt);
    %z(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_coord,'z',subs-1,subc,subt);
    
end

x=double(permute(x,[3 2 1]));
y=double(permute(y,[3 2 1]));
z=double(permute(z,[3 2 1]));


end
