function [U,t]=retrieve_seismo_sgt(snapinfo,id,M,varargin)
%
% $Date$
% $Revision$
% $LastChangedBy$

pnm_media ='../input/';
pnm_out='../output/';

%-- flags --
n=1;
while n<=nargin-3

if numel(varargin{n})==1 | ~isnumeric(varargin{n})
   switch varargin{n}
   case 'outdir'
       pnm_out=varargin{n+1}; n=n+1;
   case 'mediadir'
       pnm_media=varargin{n+1}; n=n+1;
   end
end
n=n+1;

end

disp(['outdir is ' pnm_out])
% -------------------- load data --------------------------
[rho,varnm]=gather_media(snapinfo,'rho','mediadir',pnm_media);
[lambda,varnm]=gather_media(snapinfo,'lambda','mediadir',pnm_media);
[mu,varnm]=gather_media(snapinfo,'mu','mediadir',pnm_media);

[Txx        ]=retrieve_seismo_snap(snapinfo,id,'Txx','outdir',pnm_out);
disp('  read Txx');
[Tyy        ]=retrieve_seismo_snap(snapinfo,id,'Tyy','outdir',pnm_out);
disp('  read Tyy');
[Tzz        ]=retrieve_seismo_snap(snapinfo,id,'Tzz','outdir',pnm_out);
disp('  read Tzz');
[Txy        ]=retrieve_seismo_snap(snapinfo,id,'Txy','outdir',pnm_out);
disp('  read Txy');
[Txz        ]=retrieve_seismo_snap(snapinfo,id,'Txz','outdir',pnm_out);
disp('  read Txz');
[Tyz,t      ]=retrieve_seismo_snap(snapinfo,id,'Tyz','outdir',pnm_out);
disp('  read Tyz');
Exx=zeros(size(Txx));Eyy=zeros(size(Tyy));Ezz=zeros(size(Tzz));
Exy=zeros(size(Txy));Exz=zeros(size(Txz));Eyz=zeros(size(Tyz));

E1=(lambda+mu)./( mu.*(3*lambda+2*mu) );
E2=-lambda./( 2.0*mu.*(3*lambda+2*mu) );
E3=1.0./mu;
%nt=size(Txx,ndims(Txx));
nx=size(Txx,1);ny=size(Txx,2);nz=size(Txx,3);nt=size(Txx,4);
disp(['  read (' num2str([nx ny nz nt]) ') elements'])
for k=1:nz
for j=1:ny
for i=1:nx
    Exx(i,j,k,:)=Txx(i,j,k,:).*E1(i,j,k)+Tyy(i,j,k,:).*E2(i,j,k)+Tzz(i,j,k,:).*E2(i,j,k);
    Eyy(i,j,k,:)=Txx(i,j,k,:).*E2(i,j,k)+Tyy(i,j,k,:).*E1(i,j,k)+Tzz(i,j,k,:).*E2(i,j,k);
    Ezz(i,j,k,:)=Txx(i,j,k,:).*E2(i,j,k)+Tyy(i,j,k,:).*E2(i,j,k)+Tzz(i,j,k,:).*E1(i,j,k);
    Exy(i,j,k,:)=Txy(i,j,k,:).*E3(i,j,k)*0.5;
    Exz(i,j,k,:)=Txz(i,j,k,:).*E3(i,j,k)*0.5;
    Eyz(i,j,k,:)=Tyz(i,j,k,:).*E3(i,j,k)*0.5;
end
end
end

% -------------------- represetation ------------------------
U= Exx*M(1,1)+Eyy*M(2,2)+Ezz*M(3,3) ...
  +Exy*(M(1,2)+M(2,1)) ...
  +Exz*(M(1,3)+M(3,1)) ...
  +Eyz*(M(2,3)+M(3,2)) ;
