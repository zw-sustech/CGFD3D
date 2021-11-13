function [tfr2,OrderedSurfaces]=tfrsurf(tfr,threshold,keep,trace);
% [tfr2,OrderedSurfaces]=tfrsurf(tfr,threshold,keep,trace);
% extract from a time-frequency representation the biggest energy dots
%        TFR       : time-frequency representation.
%        THRESHOLD : the energy threshold, in % 
%        KEEP      : number of dots to keep
%        TRACE     : if nonzero, the progression of the algorithm is shown
%                    (default : 0).
%
% example :
%
% N=256; 
% sig=fmlin(N,0.1,0.3)+fmlin(N,0.3,0.4)+2*fmlin(N,0.05,0.2).*amgauss(N,190,70);
% tfr=tfrwv(sig,1:N,128);
% [tfr2,OrderedSurfaces]=tfrsurf(tfr,5,3,1);
% figure(1);tfrview(tfr,sig,1:N,'tfrwv',[2 1 5 10 128 2 1 5])
% title('original tfr');
% figure(2);tfrview(tfr2,sig,1:N,'tfrwv',[2 1 5 10 128 2 1 5]);
% title('modified tfr');
% figure(3);semilogy(1:10,OrderedSurfaces(1:10),'-',1:10,OrderedSurfaces(1:10),'o');
% title('number of points of the 10 biggest dots');
%
%	See also imextract.

% F. Auger, oct 1999
%	Copyright (c) CNRS - France 1999. 
%
%	------------------- CONFIDENTIAL PROGRAM -------------------- 
%	This program can not be used without the authorization of its
%	author(s). For any comment or bug report, please send e-mail to 
%	f.auger@ieee.org 
if nargin==1, 
 threshold=5; keep=10; trace=0;
elseif nargin==2, 
 keep=10; trace=0;
elseif nargin==3, 
 trace=0; 
end;

[Nbrow,Nbcol]=size(tfr);
TheMax=max(max(tfr));
[EnergyDots,NbDots]=imextrac((tfr>=threshold*TheMax*0.01),trace);
Surfaces=zeros(1,NbDots+1);
for i=0:NbDots,
 Surfaces(i+1)=length(find(EnergyDots==i));
end;
[OrderedSurfaces,Indices]=sort(Surfaces(2:NbDots+1));
OrderedSurfaces=fliplr(OrderedSurfaces);
Indices=fliplr(Indices);

Binary=zeros(Nbrow,Nbcol);
for i=1:keep,
 DotIndice=find(EnergyDots==Indices(i));
 Binary(DotIndice)=1;
end;

tfr2=tfr.*Binary;


