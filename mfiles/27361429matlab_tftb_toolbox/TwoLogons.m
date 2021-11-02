% create a movie (or play it) showing the interferences
% between two logons

% 	F. Auger, oct 1999
%	Copyright (c) CNRS - France 1999. 
%
%	------------------- CONFIDENTIAL PROGRAM -------------------- 
%	This program can not be used without the authorization of its
%	author(s). For any comment or bug report, please send e-mail to 
%	f.auger@ieee.org 

NbRepeat=15;
if exist('twologons.mat')~=2,
 Nt=128; Nf=128; NbImages=20; 

 sig1=amgauss(Nt,Nt*3/8,40).*fmconst(Nt,0.3);
 sig2=amgauss(Nt,Nt*6/8,30).*fmconst(Nt,0.15);

 ExpJTheta=exp(j*2.0*pi*(0:NbImages-1)/NbImages);

 figure(1); 
 TheMovie=moviein(NbImages,gcf);

 for NumImage=1:NbImages,
  sig=2*sig1+ExpJTheta(NumImage)*sig2;
  tfr=tfrwv(sig,1:Nt,Nf);
  tfrview(tfr,sig,1:Nt,'tfrwv',[2 0 1 20 Nf 4 1 3]);
  Image=getframe(gcf); TheMovie(1:length(Image),NumImage)=Image;
 end;
 save twologons.mat TheMovie
 fprintf('file twologons.mat saved. So as to use it, just type\n');
 fprintf('''twologons'' in the command window.\n');
 fprintf('Don''t forget to remove it when you will no more use it.\n');
 clear TheMovie;
end;

load twologons.mat
set(gcf,'Name','Interferences between two logons');
movie(gcf,TheMovie,NbRepeat,12,[0 0 0 0])
