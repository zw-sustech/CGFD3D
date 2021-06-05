function [V,varargout] = gather_media(snapinfo,varargin)

% gather media data
% Author:   Yuanhang Huo
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.05.31

%-- flags --
nargs=nargin-1;

n=1;
while n<=nargs
    
    if numel(varargin{n})==1 | ~isnumeric(varargin{n})
        switch varargin{n}
            case {'rho','lambda','mu'}
                varnm=varargin{n};
            case 'mediadir'
                output_dir=varargin{n+1}; n=n+1;
        end
    end
    
    n=n+1;
    
end

% check path exists
if ~ exist(output_dir,'dir')
   error([mfilename ': directory ' output_dir ' does not exist']);
end

% load
mediaprefix='media';
nthd=length(snapinfo);
for n=1:nthd
    
    n_i=snapinfo(n).thisid(1); n_j=snapinfo(n).thisid(2);
    i1=snapinfo(n).indxs(1); j1=snapinfo(n).indxs(2); k1=snapinfo(n).indxs(3);
    i2=snapinfo(n).indxe(1); j2=snapinfo(n).indxe(2); k2=snapinfo(n).indxe(3);
    subs=snapinfo(n).wsubs; subc=snapinfo(n).wsubc; subt=snapinfo(n).wsubt;
    fnm_media=[output_dir,'/',mediaprefix,'_px',num2str(n_i),'_py',num2str(n_j),'.nc'];
    
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