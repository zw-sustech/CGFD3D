function [snapinfo]=locate_snap(fnm_conf,id,varargin)
% locate snap index in mpi threads
%
% $Date$
% $Revision$
% $LastChangedBy$

gsubs=[1 1 1];
gsubc=[-1 -1 -1];
gsubt=[1 1 1];
gtstart=1;gtcount=-1;gtstride=1;

%-- flags --
n=1;
while n<=nargin-2

if numel(varargin{n})==1 | ~isnumeric(varargin{n})
   switch varargin{n}
   case 'start'
        gsubs=varargin{n+1}; n=n+1;
        if length(gsubs)==4
           gtstart=gsubs(4);gsubs=gsubs(1:3);
        end
   case 'count'
        gsubc=varargin{n+1}; n=n+1;
        if length(gsubc)==4
           gtcount=gsubc(4);gsubc=gsubc(1:3);
        end
   case 'stride'
        gsubt=varargin{n+1}; n=n+1;
        if length(gsubt)==4
           gtstride=gsubt(4);gsubt=gsubt(1:3);
        end
   end
end

n=n+1;

end

% check
if ~ exist(fnm_conf,'file')
   error([mfilename ': file ' fnm_conf ' does not exist']);
end

LenFD=3;
str_cm='#';
stag=['snap_' num2str(id,'%3.3i')];

fid=fopen(fnm_conf);
conf=textscan(fid,'%s','delimiter','\n','whitespace','');
fclose(fid);
nline=size(conf{1});

for n=1:nline

str=conf{1}{n};
if isempty(str)
   continue
end
npit=findstr(str,str_cm);
if ~ isempty(npit)
   str(npit(1):end)=[];
end
str=regexprep(str,{'=','\|'},' ');
[tag,s]=strtok(str);
if isempty(tag)
   continue
end

switch tag
case 'dims'
    dims=sscanf(s,'%f',3)';
case 'ni'
    ni=sscanf(s,'%f',1);
case 'nj'
    nj=sscanf(s,'%f',1);
case 'nk'
    nk=sscanf(s,'%f',1);
case 'nt'
    nt=sscanf(s,'%f',1);
case 'GRID_ROOT'
    pnm_grid=sscanf(s,'%s',1);
case 'MEDIA_ROOT'
    pnm_media=sscanf(s,'%s',1);
case 'SOURCE_ROOT'
    pnm_src=sscanf(s,'%s',1);
case 'OUTPUT_ROOT'
    pnm_output=sscanf(s,'%s',1);
case 'number_of_recv'
    num_recv=sscanf(s,'%f',1);
case 'number_of_inline'
    num_line=sscanf(s,'%f',1);
case 'number_of_snap'
    num_snap=sscanf(s,'%f',1);
case stag
    [a,m,msg,npt]=sscanf(s,'%f',11);
    snap_subs=a(1:3)';
    snap_subc=a(4:6)';
    snap_subt=a(7:9)';
    snap_tinv=a(10);
    snap_tcnt=a(11);
    snap_stress=0;
    b=sscanf(s(npt:end),'%s');
    if ~ isempty(b)
    if strcmp(lower(b),'t')
       snap_stress=1;
    end
    end
end  % select

end

ngi=ni*dims(1);
ngj=nj*dims(2);
ngk=nk*dims(3);
ngijk=[ngi,ngj,ngk];

if id==0
    snap_subs=[  1, 1, 1 ];
    snap_subc=[ -1,-1,-1 ];
    snap_subt=[  1, 1, 1 ];
    snap_tinv=1;
    snap_tcnt=1;
end

if ~ exist('snap_subs')
    error([ 'id=',num2str(id),' doesn''t exist']);
end

%if snap_subc(1)==-1
%   snap_subc(1)=fix((ngi-snap_subs(1))/snap_subt(1))+1;
%end
%if snap_subc(2)==-1
%   snap_subc(2)=fix((ngj-snap_subs(2))/snap_subt(2))+1;
%end
%if snap_subc(3)==-1
%   snap_subc(3)=fix((ngk-snap_subs(3))/snap_subt(3))+1;
%end
%   snap_sube=snap_subs+(snap_subc-1).*snap_subt;

% reset count=-1 to actual number
indx=find(snap_subc==-1);
snap_subc(indx)=fix( (ngijk(indx)-snap_subs(indx))./snap_subt(indx) )+1;
snap_sube=snap_subs+(snap_subc-1).*snap_subt;

indx=find(gsubc==-1);
gsubc(indx)=fix( (snap_subc(indx)-gsubs(indx))./gsubt(indx) )+1;
gsube=gsubs+(gsubc-1).*gsubt;

% list of retreive points
ilist=snap_subs(1):snap_subt(1):snap_sube(1);
jlist=snap_subs(2):snap_subt(2):snap_sube(2);
klist=snap_subs(3):snap_subt(3):snap_sube(3);

rlist=gsubs(1):gsubt(1):gsube(1);
slist=gsubs(2):gsubt(2):gsube(2);
tlist=gsubs(3):gsubt(3):gsube(3);

ulist=ilist(rlist);
vlist=jlist(slist);
wlist=klist(tlist);

% find the mpi threads scale
n_i1=fix( (ulist(1)  -1)/ni);
n_i2=fix( (ulist(end)-1)/ni);
n_j1=fix( (vlist(1)  -1)/nj);
n_j2=fix( (vlist(end)-1)/nj);
n_k1=fix( (wlist(1)  -1)/nk);
n_k2=fix( (wlist(end)-1)/nk);

% locate
nthd=0;
for n_i=n_i1:n_i2
for n_j=n_j1:n_j2
for n_k=n_k1:n_k2
    nthd=nthd+1;
    snapinfo(nthd).thisid=[n_i,n_j,n_k];
    i=find(ulist>=ni*n_i+1 & ulist<=ni*(n_i+1));
    j=find(vlist>=nj*n_j+1 & vlist<=nj*(n_j+1));
    k=find(wlist>=nk*n_k+1 & wlist<=nk*(n_k+1));
    snapinfo(nthd).indxs=[i(1),j(1),k(1)];
    snapinfo(nthd).indxe=[i(end),j(end),k(end)];
    snapinfo(nthd).indxc=snapinfo(nthd).indxe-snapinfo(nthd).indxs+1;

    snapinfo(nthd).wsubs=[ ulist(i(1))-ni*n_i, ...
                           vlist(j(1))-nj*n_j, ...
                           wlist(k(1))-nk*n_k]+LenFD;
    snapinfo(nthd).wsubc=snapinfo(nthd).indxc;
    snapinfo(nthd).wsubt=snap_subt.*gsubt;

    r=find(ilist>=ni*n_i+1 & ilist<=ni*(n_i+1));
    s=find(jlist>=nj*n_j+1 & jlist<=nj*(n_j+1));
    t=find(klist>=nk*n_k+1 & klist<=nk*(n_k+1));
    snapinfo(nthd).subs=[ (ulist(i(1))-ilist(r(1)))/snap_subt(1)+1, ...
                          (vlist(j(1))-jlist(s(1)))/snap_subt(2)+1, ...
                          (wlist(k(1))-klist(t(1)))/snap_subt(3)+1 ];
    snapinfo(nthd).subc=snapinfo(nthd).indxc;
    snapinfo(nthd).subt=gsubt;

    snapinfo(nthd).tinv=snap_tinv;
    snapinfo(nthd).tcnt=snap_tcnt;

    snapinfo(nthd).ttriple=[gtstart,gtcount,gtstride];
end
end
end

