function [V,varargout]=gather_image(snapinfo,id,nlayer,varargin)
% gather snap data
%
% $Date: 2009-01-07 22:26:46 -0500 (Wed, 07 Jan 2009) $
% $Revision: 503 $
% $LastChangedBy: zhangw $

varnm = 'image';
fnm_prefix='image_';

pnm_out='../RTMoutput/';

%-- flags --
nargs=nargin-3
n=1;

while n<=nargs

if numel(varargin{n})==1 | ~isnumeric(varargin{n})
   switch varargin{n}
   case {'Vx','Vy','Vz'}
       varnm=varargin{n};
       fnm_prefix='vel_';
   case {'Ax','Ay','Az'}
       varnm=varargin{n};
       fnm_prefix='acce_';
   case {'Dx','Dy','Dz'}
       varnm=varargin{n};
       fnm_prefix='disp_';
   case {'Txx','Tyy','Tzz','Txy','Txz','Tyz'}
       varnm=varargin{n};
       fnm_prefix='sgt_';
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

nofnc=fix( (nlayer-1)/snapinfo(1).tcnt )+1;
ninnc=mod( (nlayer-1),snapinfo(1).tcnt );
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
    tdim=nc_getdiminfo(fnm_snap,'time');
    if tdim.Length==0 | ninnc-1>=tdim.Length
       error([num2str(ninnc+1) 'th layer is beyond current time dim(' ...
            num2str(tdim.Length) ') in ' fnm_snap]);
       %return
    end
    %disp([num2str(ninnc) 'th layer in ',fnm_snap]);
    %-- get data
    V(k1:k2,j1:j2,i1:i2)=nc_varget(fnm_snap,varnm, ...
          [ninnc,fliplr(subs)-1],[1,fliplr(subc)],[1,fliplr(subt)]);
    if nargout>=2 & n==1, t=nc_varget(fnm_snap,'time',[ninnc],[1]); end
end

V=permute(V,[3,2,1]);

if nargout>=2, varargout(1)={t}; end
if nargout>=3,varargout(2)={varnm}; end
