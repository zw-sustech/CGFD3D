function [V,varargout]=retrieve_seismo_snap(snapinfo,id,varargin)
% retreive seismogram from snap nc files
%
% $Date$
% $Revision$
% $LastChangedBy$

varnm = 'Vz';
fnm_prefix='snap_';
pnm_out='../output/';

nout=nargout-1;

%-- flags --
n=1;
while n<=nargin-2

if numel(varargin{n})==1 | ~isnumeric(varargin{n})
   switch varargin{n}
   case 'outdir'
       pnm_out=varargin{n+1}; n=n+1;
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
   end
end
n=n+1;

end

% check
if ~ exist(pnm_out,'dir')
   error([mfilename ': directory ' pnm_out ' does not exist']);
end

nthd=length(snapinfo);

for n=1:nthd
    n_i=snapinfo(n).thisid(1);n_j=snapinfo(n).thisid(2);n_k=snapinfo(n).thisid(3);
    i1=snapinfo(n).indxs(1);j1=snapinfo(n).indxs(2);k1=snapinfo(n).indxs(3);
    i2=snapinfo(n).indxe(1);j2=snapinfo(n).indxe(2);k2=snapinfo(n).indxe(3);
    subs=snapinfo(n).subs;
    subc=snapinfo(n).subc;
    subt=snapinfo(n).subt;

nofnc=0; nt=0;
ninnc=snapinfo(n).ttriple(1);
ncount=snapinfo(n).ttriple(2);
nstride=snapinfo(n).ttriple(3);
while 1
    nofnc=nofnc+1;
    fnm_snap=get_fnm_snap(pnm_out,fnm_prefix,n_i,n_j,n_k,id,nofnc);
    if ~ exist(fnm_snap)
       break
    end

    tdim=nc_getdiminfo(fnm_snap,'time');
    if tdim.Length==0
       break
    end
    LenT=tdim.Length;
    ntinc=fix( (LenT-ninnc)/nstride )+1;

    if ncount ~= -1 & ntinc+nt>=ncount
       ntinc=ncount-nt;
    end
    %-- get data
    if nout>0
        t(nt+1:nt+ntinc)=nc_varget(fnm_snap,'time' ,[ninnc-1],[ntinc],[nstride]);
    end
    V(nt+1:nt+ntinc,k1:k2,j1:j2,i1:i2)=nc_varget(fnm_snap,varnm, ...
          [ninnc-1,fliplr(subs)-1],[ntinc,fliplr(subc)],[nstride,fliplr(subt)]);

    nt=nt+ntinc;
    ninnc=nstride-(LenT-(ninnc+(ntinc-1)*nstride));
    if ncount ~= -1 & nt>=ncount
       break
    end
end

end

V=permute(V,[4,3,2,1]);
if nout>0, varargout{1}=t; end
if nout>1, varargout{2}=varnm; end
