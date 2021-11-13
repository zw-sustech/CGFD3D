function [px,py,pz]=fun_border4slice(xd,yd,zd,scl_xlim,scl_ylim,MDLA,MDLO,varargin)
%
% fun_border4slice: Find the border lines coordinate.
%
% Usage: [px,py,pz]=fun_border4slice(xd,yd,zd,scl_xlim,scl_ylim,MDLA,MDLO,varargin)
%

% Major ChangeLog:
%   2009-01-09 Wei Zhang
%     * Added help information, but uncomplete.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date$
% $Revision$
% $LastChangedBy$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_hide=0;

% ----------------------- parameter -----------------------
args=varargin; nargs=numel(args);

n=1;
while n<=nargs 

if ~ischar(args{n}), n=n+1; continue; end

switch args{n}
case 'hide'
   flag_hide=1;
end

n=n+1;
end

% ----------------------- border -----------------------
nslc=numel(xd);
ptx=[scl_xlim(1),scl_xlim(2),scl_xlim(2),scl_xlim(1),scl_xlim(1)];
pty=[scl_ylim(1),scl_ylim(1),scl_ylim(2),scl_ylim(2),scl_ylim(1)];

% -- vertical plane --
for n=1:2
px{n}=[xd{n}(:,1)',xd{n}(end,:),xd{n}(end:-1:1,end)',xd{n}(1,end:-1:1)];
py{n}=[yd{n}(:,1)',yd{n}(end,:),yd{n}(end:-1:1,end)',yd{n}(1,end:-1:1)];
pz{n}=[zd{n}(:,1)',zd{n}(end,:),zd{n}(end:-1:1,end)',zd{n}(1,end:-1:1)];

indx11=find(px{n}<scl_xlim(1)); indx12=find(px{n}>scl_xlim(2));
indx21=find(py{n}<scl_ylim(1)); indx22=find(py{n}>scl_ylim(2));
px{n}([indx11,indx12,indx21,indx22])=[];
py{n}([indx11,indx12,indx21,indx22])=[];
pz{n}([indx11,indx12,indx21,indx22])=[];
if px{n}(end)~=px{n}(1) | py{n}(end)~=py{n}(1) | pz{n}(end)~=pz{n}(1)
   px{n}(end+1)=px{n}(1);py{n}(end+1)=py{n}(1);pz{n}(end+1)=pz{n}(1);
end
end

% -- horizontal plane --
for n=3:nslc

px{n}=[xd{n}(1,1),xd{n}(end,1),xd{n}(end,end),xd{n}(1,end),xd{n}(1,1)];
py{n}=[yd{n}(1,1),yd{n}(end,1),yd{n}(end,end),yd{n}(1,end),yd{n}(1,1)];
pz{n}=[zd{n}(1,1),zd{n}(end,1),zd{n}(end,end),zd{n}(1,end),zd{n}(1,1)];

M=inpolygon(ptx(1:4),pty(1:4),px{n},py{n}); indx=find(M);
[xint,yint]=polyxpoly(ptx,pty,px{n},py{n});
if length(indx)==1
   px{n}(2:4)=[xint(1),ptx(indx),xint(2)];
   py{n}(2:4)=[yint(1),pty(indx),yint(2)];
elseif length(indx)==2
   px{n}(6)=px{n}(5);py{n}(6)=py{n}(5);pz{n}(6)=pz{n}(5);
   px{n}(2:5)=[xint(1),ptx(indx(1)),ptx(indx(2)),xint(2)];
   py{n}(2:5)=[yint(1),pty(indx(1)),pty(indx(2)),yint(2)];
   if polyarea(px{n},py{n}) ...
     <polyarea([px{n}(1:2),px{n}(4:-1:3),px{n}(5:6)],[py{n}(1:2),py{n}(4:-1:3),py{n}(5:6)])
     px{n}(3:4)=px{n}(4:-1:3);py{n}(3:4)=py{n}(4:-1:3);
   end
elseif length(indx)>=3
   error('why morn than two points inpolygon a othogonal poly plane?');
else
   px{n}=[px{n}(1),xint(1),xint(2),px{n}(end)];
   py{n}=[py{n}(1),yint(1),yint(2),py{n}(end)];
   pz{n}(end)=[];
end

end

% -- hide horizontal plane --
if flag_hide
   for n=4:min(5,nslc)
       if py{n}(end-1)>py{n}(2)
          py{n}=py{n}(end:-1:1); px{n}=px{n}(end:-1:1);
       end
   py{n}(end)=[];px{n}(end)=[];pz{n}(end)=[];
   end
end

% center line
px{end+1}=[MDLO MDLO];
py{end+1}=[MDLA MDLA];
pz{end+1}=[zd{1}(1,1),zd{1}(end,end)];

