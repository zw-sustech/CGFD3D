function [V,varargout] = gather_media(mediainfo,varnm,media_dir)

% gather media data
% Author:   Yuanhang Huo
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.05.31

%nargs=nargin-1;

%n=1;
%while n<=nargs
%    
%    if numel(varargin{n})==1 | ~isnumeric(varargin{n})
%        switch varargin{n}
%            case {'rho','lambda','mu','kappa'}
%                varnm=varargin{n};
%            case 'mediadir'
%                media_dir=varargin{n+1}; n=n+1;
%        end
%    end
%    
%    n=n+1;
%    
%end

% check path exists
if ~ exist(media_dir,'dir')
   error([mfilename ': directory ' media_dir ' does not exist']);
end

% load
mediaprefix='media';
nthd=length(mediainfo);
for n=1:nthd
    
    n_i=mediainfo(n).thisid(1); n_j=mediainfo(n).thisid(2);
    i1=mediainfo(n).indxs(1); j1=mediainfo(n).indxs(2); k1=mediainfo(n).indxs(3);
    i2=mediainfo(n).indxe(1); j2=mediainfo(n).indxe(2); k2=mediainfo(n).indxe(3);
    subs=mediainfo(n).wsubs; subc=mediainfo(n).wsubc; subt=mediainfo(n).wsubt;
    fnm_media=[media_dir,'/',mediaprefix,'_px',num2str(n_i),'_py',num2str(n_j),'.nc'];
    
    if ~ exist(fnm_media,'file')
       error([mfilename ': file ' fnm_media 'does not exist']);
    end

    V(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_media,varnm, ...
                            fliplr(subs)-1,fliplr(subc),fliplr(subt));
                               
end

V=double(permute(V,[3,2,1]));

% pack varargout
if nargout>=2
    varargout(1)={varnm};
end


end
