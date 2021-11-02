% create a movie (or play it) showing the trade-off
% for the pweudo Wigner-Ville Distribution

% 	F. Auger, oct 1999
%	Copyright (c) CNRS - France 1999. 
%
%	------------------- CONFIDENTIAL PROGRAM -------------------- 
%	This program can not be used without the authorization of its
%	author(s). For any comment or bug report, please send e-mail to 
%	f.auger@ieee.org 

NbRepeat=6;
if exist('pwvmovie.mat')~=2,
 Nt=256; Nf=256; NbImages=30; 

 sig1=amgauss(Nt,Nt*2/32, 10).*fmconst(Nt,0.20,Nt*1/16); sig1=sig1/norm(sig1); 
 sig2=amgauss(Nt,Nt*5/32, 10).*fmconst(Nt,0.20,Nt*3/16); sig2=sig2/norm(sig2); 
 sig3=amgauss(Nt,Nt*9/16,110).*fmconst(Nt,0.35,Nt*9/16); sig3=sig3/norm(sig3); 
 sig4=amgauss(Nt,Nt*9/16,110).*fmconst(Nt,0.40,Nt*9/16); sig4=sig4/norm(sig4); 
 sig=sig1+sig2+2*(sig3+sig4);

 Lh=round(linspace(0,1,NbImages).^(2)*126)+1;

 figure(1); 
 TheMovie=moviein(NbImages,gcf);

 for NumImage=1:NbImages,
  h=window(2*Lh(NumImage)+1,'hanning');
  tfr=tfrpwv(sig,1:Nt,Nf,h);
  tfrview(tfr,sig,1:Nt,'TFRPWV',[2 0 2 20 Nf 4 1 3],h);
  Image=getframe(gcf); TheMovie(1:length(Image),NumImage)=Image;
 end;
 save pwvmovie.mat TheMovie
 fprintf('file pwvmovie.mat saved. So as to use it, just type\n');
 fprintf('''pwvmovie'' in the command window.\n');
 fprintf('Don''t forget to remove it when you will no more use it.\n');
 clear TheMovie;
end;

load pwvmovie.mat
set(gcf,'Name','Interferences between 3 signal components');
movie(gcf,TheMovie,-NbRepeat,6,[0 0 0 0])
