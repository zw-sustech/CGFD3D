function [v,varargout]=retrieve_seismo(seismoinfo,varnm,varargin)
% retreive seismogram on individual receivers.
% id=0 means recv type, otherwise line type
%
% $Date$
% $Revision$
% $LastChangedBy$

pnm_out='../output/';

%-- flags --
narg = max(nargin,1)-2; n=1;
while n<=narg

if numel(varargin{n})==1 | ~isnumeric(varargin{n})
   switch varargin{n}
   case 'outdir'
       pnm_out=varargin{n+1}; n=n+1;
   end
end
n=n+1;

end

% check
if ~ exist(pnm_out,'dir')
   error([mfilename ': directory ' pnm_out ' does not exist']);
end

fnm_nc=get_fnm_seismo(pnm_out,seismoinfo.n_i,seismoinfo.n_j,seismoinfo.n_k);

v=nc_varget(fnm_nc,varnm,[0 seismoinfo.npt-1],[-1,1]);
if nargout>1
   varargout(1)={nc_varget(fnm_nc,'time')};
end
