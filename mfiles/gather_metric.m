function [V,varargout] = gather_metric(metricinfo,varargin)

% gather metric data
% Author:   Yuanhang Huo
% Email:    yhhuo@mail.ustc.edu.cn
% Date:     2021.06.06

nargs=nargin-1;

n=1;
while n<=nargs
    
    if numel(varargin{n})==1 | ~isnumeric(varargin{n})
        switch varargin{n}
            case {'jac','xi_x','xi_y','xi_z',...
                  'eta_x','eta_y','eta_z',...
                  'zeta_x','zeta_y','zeta_z'}
                varnm=varargin{n};
            case 'metricdir'
                metric_dir=varargin{n+1}; n=n+1;
        end
    end
    
    n=n+1;
    
end

% check path exists
if ~ exist(metric_dir,'dir')
   error([mfilename ': directory ' metric_dir ' does not exist']);
end

% load
metricprefix='metric';
nthd=length(metricinfo);
for n=1:nthd
    
    n_i=metricinfo(n).thisid(1); n_j=metricinfo(n).thisid(2);
    i1=metricinfo(n).indxs(1); j1=metricinfo(n).indxs(2); k1=metricinfo(n).indxs(3);
    i2=metricinfo(n).indxe(1); j2=metricinfo(n).indxe(2); k2=metricinfo(n).indxe(3);
    subs=metricinfo(n).wsubs; subc=metricinfo(n).wsubc; subt=metricinfo(n).wsubt;
    fnm_metric=[metric_dir,'/',metricprefix,'_px',num2str(n_i),'_py',num2str(n_j),'.nc'];
    
    if ~ exist(fnm_metric,'file')
       error([mfilename ': file ' fnm_metric 'does not exist']);
    end

    V(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_metric,varnm, ...
                            fliplr(subs)-1,fliplr(subc),fliplr(subt));
                               
end

V=double(permute(V,[3,2,1]));

% pack varargout
if nargout>=2
    varargout(1)={varnm};
end


end