function flag=istfr1(method);
% flag=istfr1(method) returns true is method is a
% time frequency representation of type 1 (positive and negative frequencies).
%	See also istfr2, istfraff.

%	F. Auger, may 98
%	Copyright (c) CNRS - France 1998. 
%
%	------------------- CONFIDENTIAL PROGRAM -------------------- 
%	This program can not be used without the authorization of its
%	author(s). For any comment or bug report, please send e-mail to 
%	f.auger@ieee.org 

method=upper(method);
if strcmp(method,'TFRPMH'  )| strcmp(method,'TFRRPMH' )| ...
   strcmp(method,'TFRSP'   )| strcmp(method,'TFRRSP'  )| ...
   strcmp(method,'TFRPPAGE')| strcmp(method,'TFRRPPAG')| ...
   strcmp(method,'TFRMHS'  )| strcmp(method,'TFRRGAB' )| ...
   strcmp(method,'TFRMH'   )| strcmp(method,'TFRMMCE' )| ...
   strcmp(method,'TFRRMSC' )| strcmp(method,'TFRPAGE' )| ...
   strcmp(method,'TFRGABOR')| strcmp(method,'TFRRI'   )| ...
   strcmp(method,'TFRMSC'  )| strcmp(method,'TYPE1'   )| ...
   strcmp(method,'TFRSTFT' ),
 flag=1;
else
 flag=0;
end;

