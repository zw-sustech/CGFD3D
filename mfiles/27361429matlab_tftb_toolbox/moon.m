% create a movie (or play it) showing the interferences
% between a sinusoid, a logon and a rotating logon

% 	F. Auger, oct 1999
%	Copyright (c) CNRS - France 1999. 
%
%	------------------- CONFIDENTIAL PROGRAM -------------------- 
%	This program can not be used without the authorization of its
%	author(s). For any comment or bug report, please send e-mail to 
%	f.auger@ieee.org 

NbRepeat=15;
if exist('moon.mat')~=2,
 Nt=128; Nf=128; NbImages=30; 

 sig1=fmconst(Nt,0.03);
 sig2=amgauss(Nt,Nt/2,30).*fmconst(Nt,0.3);

 figure(1); 
 TheMovie=moviein(NbImages,gcf);

 Theta=2.0*pi*(0:NbImages-1)/NbImages;

 for NumImage=1:NbImages,
  sig3=amgauss(Nt,Nt*(0.5+0.3*cos(Theta(NumImage))),30).*fmconst(Nt,0.3+0.15*sin(Theta(NumImage)));
  sig=0.5*sig1+0.5*sig2+sig3;
  tfr=tfrwv(sig,1:Nt,Nf);
  tfrview(tfr,sig,1:Nt,'tfrwv',[2 0 1 20 Nf 4 1 3 0 0.5]);
  Image=getframe(gcf); TheMovie(1:length(Image),NumImage)=Image;
 end;
 save moon.mat TheMovie
 fprintf('file moon.mat saved. So as to use it, just type\n');
 fprintf('''moon'' in the command window.\n');
 fprintf('Don''t forget to remove it when you will no more use it.\n');
 clear TheMovie;
end;

load moon.mat
set(gcf,'Name','Interferences between 3 signal components');
movie(gcf,TheMovie,NbRepeat,12,[0 0 0 0])
