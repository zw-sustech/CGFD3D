function [V,varargout]=gather_dist(snapinfo,id,varargin)
% gather snap data
%
% $Date$
% $Revision$
% $LastChangedBy$

varnm = 'Vz';
fnm_prefix='Vmax_';

pnm_out='../output/';

%-- flags --
nargs=nargin-2;

n=1;
while n<=nargs

if numel(varargin{n})==1 | ~isnumeric(varargin{n})
   switch varargin{n}
   case {'Vx','Vy','Vz','Vn','Ve','Vh','Va'}
       varnm=varargin{n};
       fnm_prefix='Vmax_';
   case {'Ax','Ay','Az','An','Ae','Ah','Aa'}
       varnm=varargin{n};
       fnm_prefix='Amax_';
   case {'Dx','Dy','Dz','Dn','De','Dh','Da'}
       varnm=varargin{n};
       fnm_prefix='Dmax_';
   case {'Kapx','Kaqx','Kbpx','Kbqx', ...
         'Kapy','Kaqy','Kbpy','Kbqy', ...
         'Kapz','Kaqz','Kbpz','Kbqz' }
       varnm=varargin{n};
       fnm_prefix=['kernel_V' varnm(end) '_'];
       varnm=varnm(1:end-1);
   case {'kernel_phase_Vp','kernel_amplitude_Vp', ...
         'kernel_phase_Vs','kernel_amplitude_Vs'}
       varnm=varargin{n};
       fnm_prefix=['kernel_'];
   case {'phase_Vp','amplitude_Vp', ...
         'phase_Vs','amplitude_Vs'}
       varnm=varargin{n};
       fnm_prefix=['kernel_'];
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

nthd=length(snapinfo);

nofnc=0;
for n=1:nthd
    n_i=snapinfo(n).thisid(1);n_j=snapinfo(n).thisid(2);n_k=snapinfo(n).thisid(3);
    i1=snapinfo(n).indxs(1);j1=snapinfo(n).indxs(2);k1=snapinfo(n).indxs(3);
    i2=snapinfo(n).indxe(1);j2=snapinfo(n).indxe(2);k2=snapinfo(n).indxe(3);
    subs=snapinfo(n).subs;
    subc=snapinfo(n).subc;
    subt=snapinfo(n).subt;
    fnm_snap=get_fnm_snap(pnm_out,fnm_prefix,n_i,n_j,n_k,id,nofnc);
    if ~ exist(fnm_snap)
       error([mfilename ': file ',fnm_snap, ' does not exist']);
    end
    %-- get data
    %V(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_snap,varnm,[0,0,0],[subc]);
    V(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_snap,varnm, ...
          fliplr(subs)-1,fliplr(subc),fliplr(subt));
end

V=permute(V,[3,2,1]);

% pack varargout
if nargout>=2, varargout(1)={varnm}; end

