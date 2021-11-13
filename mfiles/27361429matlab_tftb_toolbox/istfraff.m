function flag=istfraff(method);
% flag=istfr2(method) returns true is method is an affine
% time frequency representation.
%	See also istfr1, istfr2.

%	F. Auger, may 98
%	Copyright (c) CNRS - France 1998. 
%
%	------------------- CONFIDENTIAL PROGRAM -------------------- 
%	This program can not be used without the authorization of its
%	author(s). For any comment or bug report, please send e-mail to 
%	f.auger@ieee.org 

method=upper(method);
if strcmp(method,'TFRASPW' ) | strcmp(method,'TFRSCALO') | ...
   strcmp(method,'TFRDFLA' ) | strcmp(method,'TFRSPAW' ) | ...
   strcmp(method,'TFRUNTER') | strcmp(method,'TFRBERT' ) | ...
   strcmp(method,'TFRSPBK' ),
 flag=1;
else
 flag=0;
end;
