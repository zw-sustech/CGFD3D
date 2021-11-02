%LH    list SAC header
%
%    Read or set matlab variables to SAC header variables from
%    SAC files read in to matlab with rsac.m
%    
%    Examples:
%
%    To list all defined header variables in the file KATH:
%    lh(KATH)  
%
%    To assign the SAC variable DELTA from station KATH to
%    the matlab variable dt:
%
%    dt = lh(KATH,'DELTA');
%
%    To assign the SAC variables STLA and STLO from station KATH
%    to the matlab variables lat and lon:
%
%    [lat,lon] = lh(KATH,'STLA','STLO')
%
%    by Michael Thorne (4/2004)  mthorne@asu.edu
%
%    See also:  RSAC, CH, BSAC, WSAC 

function [varargout] = lh(file,varargin);

% first test to see if the file is indeed a sacfile
%---------------------------------------------------------------------------
if (file(303,3)~=77 & file(304,3)~=73 & file(305,3)~=75 & file(306,3)~=69)
  error('Specified Variable is not in SAC format ...')
end

h(1:306) = file(1:306,3); 


% read real header variables
%---------------------------------------------------------------------------
DELTA = h(1);
if (h(1) ~= -12345 & nargin == 1); disp(sprintf('DELTA      = %0.8g',h(1))); end
DEPMIN = h(2);
if (h(2) ~= -12345 & nargin == 1); disp(sprintf('DEPMIN     = %0.8g',h(2))); end
DEPMAX = h(3);
if (h(3) ~= -12345 & nargin == 1); disp(sprintf('DEPMAX     = %0.8g',h(3))); end
SCALE = h(4);
if (h(4) ~= -12345 & nargin == 1);  disp(sprintf('SCALE      = %0.8g',h(4))); end
ODELTA = h(5);
if (h(5) ~= -12345 & nargin == 1);  disp(sprintf('ODELTA     = %0.8g',h(5))); end
B = h(6);
if (h(6) ~= -12345 & nargin == 1); disp(sprintf('B          = %0.8g',h(6))); end
E = h(7);
if (h(7) ~= -12345 & nargin == 1); disp(sprintf('E          = %0.8g',h(7))); end
O = h(8);
if (h(8) ~= -12345 & nargin == 1); disp(sprintf('O          = %0.8g',h(8))); end 
A = h(9);
if (h(9) ~= -12345 & nargin == 1); disp(sprintf('A          = %0.8g',h(9))); end 
T0 = h(11);
if (h(11) ~= -12345 & nargin == 1); disp(sprintf('T0         = %0.8g',h(11))); end
T1 = h(12);
if (h(12) ~= -12345 & nargin == 1); disp(sprintf('T1         = %0.8g',h(12))); end
T2 = h(13);
if (h(13) ~= -12345 & nargin == 1); disp(sprintf('T2         = %0.8g',h(13))); end
T3 = h(14);
if (h(14) ~= -12345 & nargin == 1); disp(sprintf('T3         = %0.8g',h(14))); end
T4 = h(15);
if (h(15) ~= -12345 & nargin == 1); disp(sprintf('T4         = %0.8g',h(15))); end
T5 = h(16);
if (h(16) ~= -12345 & nargin == 1); disp(sprintf('T5         = %0.8g',h(16))); end
T6 = h(17);
if (h(17) ~= -12345 & nargin == 1); disp(sprintf('T6         = %0.8g',h(17))); end
T7 = h(18);
if (h(18) ~= -12345 & nargin == 1); disp(sprintf('T7         = %0.8g',h(18))); end
T8 = h(19);
if (h(19) ~= -12345 & nargin == 1); disp(sprintf('T8         = %0.8g',h(19))); end
T9 = h(20);
if (h(20) ~= -12345 & nargin == 1); disp(sprintf('T9         = %0.8g',h(20))); end
F = h(21);
if (h(21) ~= -12345 & nargin == 1); disp(sprintf('F          = %0.8g',h(21))); end 
RESP0 = h(22);
if (h(22) ~= -12345 & nargin == 1); disp(sprintf('RESP0      = %0.8g',h(22))); end
RESP1 = h(23);
if (h(23) ~= -12345 & nargin == 1); disp(sprintf('RESP1      = %0.8g',h(23))); end
RESP2 = h(24);
if (h(24) ~= -12345 & nargin == 1); disp(sprintf('RESP2      = %0.8g',h(24))); end
RESP3 = h(25);
if (h(25) ~= -12345 & nargin == 1); disp(sprintf('RESP3      = %0.8g',h(25))); end
RESP4 = h(26);
if (h(26) ~= -12345 & nargin == 1); disp(sprintf('RESP4      = %0.8g',h(26))); end
RESP5 = h(27);
if (h(27) ~= -12345 & nargin == 1); disp(sprintf('RESP5      = %0.8g',h(27))); end
RESP6 = h(28);
if (h(28) ~= -12345 & nargin == 1); disp(sprintf('RESP6      = %0.8g',h(28))); end
RESP7 = h(29);
if (h(29) ~= -12345 & nargin == 1); disp(sprintf('RESP7      = %0.8g',h(29))); end
RESP8 = h(30);
if (h(30) ~= -12345 & nargin == 1); disp(sprintf('RESP8      = %0.8g',h(30))); end
RESP9 = h(31);
if (h(31) ~= -12345 & nargin == 1); disp(sprintf('RESP9      = %0.8g',h(31))); end
STLA = h(32);
if (h(32) ~= -12345 & nargin == 1); disp(sprintf('STLA       = %0.8g',h(32))); end
STLO = h(33);
if (h(33) ~= -12345 & nargin == 1); disp(sprintf('STLO       = %0.8g',h(33))); end
STEL = h(34);
if (h(34) ~= -12345 & nargin == 1); disp(sprintf('STEL       = %0.8g',h(34))); end
STDP = h(35);
if (h(35) ~= -12345 & nargin == 1); disp(sprintf('STDP       = %0.8g',h(35))); end
EVLA = h(36);
if (h(36) ~= -12345 & nargin == 1); disp(sprintf('EVLA       = %0.8g',h(36))); end
EVLO = h(37);
if (h(37) ~= -12345 & nargin == 1); disp(sprintf('EVLO       = %0.8g',h(37))); end
EVEL = h(38);
if (h(38) ~= -12345 & nargin == 1); disp(sprintf('EVEL       = %0.8g',h(38))); end
EVDP = h(39);
if (h(39) ~= -12345 & nargin == 1); disp(sprintf('EVDP       = %0.8g',h(39))); end
MAG = h(40);
if (h(40) ~= -12345 & nargin == 1); disp(sprintf('MAG        = %0.8g',h(40))); end 
USER0 = h(41);
if (h(41) ~= -12345 & nargin == 1); disp(sprintf('USER0      = %0.8g',h(41))); end
USER1 = h(42);
if (h(42) ~= -12345 & nargin == 1); disp(sprintf('USER1      = %0.8g',h(42))); end
USER2 = h(43);
if (h(43) ~= -12345 & nargin == 1); disp(sprintf('USER2      = %0.8g',h(43))); end
USER3 = h(44);
if (h(44) ~= -12345 & nargin == 1); disp(sprintf('USER3      = %0.8g',h(44))); end
USER4 = h(45);
if (h(45) ~= -12345 & nargin == 1); disp(sprintf('USER4      = %0.8g',h(45))); end
USER5 = h(46);
if (h(46) ~= -12345 & nargin == 1); disp(sprintf('USER5      = %0.8g',h(46))); end
USER6 = h(47);
if (h(47) ~= -12345 & nargin == 1); disp(sprintf('USER6      = %0.8g',h(47))); end
USER7 = h(48);
if (h(48) ~= -12345 & nargin == 1); disp(sprintf('USER7      = %0.8g',h(48))); end
USER8 = h(49);
if (h(49) ~= -12345 & nargin == 1);  disp(sprintf('USER8     = %0.8g',h(49))); end
USER9 = h(50);
if (h(50) ~= -12345 & nargin == 1); disp(sprintf('USER9      = %0.8g',h(50))); end
DIST = h(51);
if (h(51) ~= -12345 & nargin == 1); disp(sprintf('DIST       = %0.8g',h(51))); end 
AZ = h(52);
if (h(52) ~= -12345 & nargin == 1); disp(sprintf('AZ         = %0.8g',h(52))); end 
BAZ = h(53);
if (h(53) ~= -12345 & nargin == 1); disp(sprintf('BAZ        = %0.8g',h(53))); end
GCARC = h(54);
if (h(54) ~= -12345 & nargin == 1); disp(sprintf('GCARC      = %0.8g',h(54))); end
DEPMEN = h(57);
if (h(57) ~= -12345 & nargin == 1); disp(sprintf('DEPMEN     = %0.8g',h(57))); end
CMPAZ = h(58);
if (h(58) ~= -12345 & nargin == 1); disp(sprintf('CMPAZ      = %0.8g',h(58))); end
CMPINC = h(59);
if (h(59) ~= -12345 & nargin == 1); disp(sprintf('CMPINC     = %0.8g',h(59))); end
XMINIMUM = h(60);
if (h(60) ~= -12345 & nargin == 1); disp(sprintf('XMINIMUM   = %0.8g',h(60))); end
XMAXIMUM = h(61);
if (h(61) ~= -12345 & nargin == 1); disp(sprintf('XMAXIMUM   = %0.8g',h(61))); end
YMINIMUM = h(62);
if (h(62) ~= -12345 & nargin == 1); disp(sprintf('YMINIMUM   = %0.8g',h(62))); end
YMAXIMUM = h(63);
if (h(63) ~= -12345 & nargin == 1); disp(sprintf('YMAXIMUM   = %0.8g',h(63))); end

% read integer header variables
%---------------------------------------------------------------------------
NZYEAR = round(h(71));
if (h(71) ~= -12345 & nargin == 1); disp(sprintf('NZYEAR     = %d',h(71))); end
NZJDAY = round(h(72));
if (h(72) ~= -12345 & nargin == 1); disp(sprintf('NZJDAY     = %d',h(72))); end
NZHOUR = round(h(73));
if (h(73) ~= -12345 & nargin == 1); disp(sprintf('NZHOUR     = %d',h(73))); end
NZMIN = round(h(74));
if (h(74) ~= -12345 & nargin == 1); disp(sprintf('NZMIN      = %d',h(74))); end
NZSEC = round(h(75));
if (h(75) ~= -12345 & nargin == 1); disp(sprintf('NZSEC      = %d',h(75))); end
NZMSEC = round(h(76));
if (h(76) ~= -12345 & nargin == 1); disp(sprintf('NZMSEC     = %d',h(76))); end
NVHDR = round(h(77));
if (h(77) ~= -12345 & nargin == 1); disp(sprintf('NVHDR      = %d',h(77))); end
NORID = round(h(78));
if (h(78) ~= -12345 & nargin == 1); disp(sprintf('NORID      = %d',h(78))); end
NEVID = round(h(79));
if (h(79) ~= -12345 & nargin == 1); disp(sprintf('NEVID      = %d',h(79))); end
NPTS = round(h(80));
if (h(80) ~= -12345 & nargin == 1); disp(sprintf('NPTS       = %d',h(80))); end
NWFID = round(h(82));
if (h(82) ~= -12345 & nargin == 1); disp(sprintf('NWFID      = %d',h(82))); end
NXSIZE = round(h(83));
if (h(83) ~= -12345 & nargin == 1); disp(sprintf('NXSIZE     = %d',h(83))); end
NYSIZE = round(h(84));
if (h(84) ~= -12345 & nargin == 1); disp(sprintf('NYSIZE     = %d',h(84))); end
IFTYPE = round(h(86));
if (h(86) ~= -12345 & nargin == 1); disp(sprintf('IFTYPE     = %d',h(86))); end
IDEP = round(h(87));
if (h(87) ~= -12345 & nargin == 1); disp(sprintf('IDEP       = %d',h(87))); end
IZTYPE = round(h(88));
if (h(88) ~= -12345 & nargin == 1); disp(sprintf('IZTYPE     = %d',h(88))); end
IINST = round(h(90));
if (h(90) ~= -12345 & nargin == 1); disp(sprintf('IINST      = %d',h(90))); end
ISTREG = round(h(91));
if (h(91) ~= -12345 & nargin == 1); disp(sprintf('ISTREG     = %d',h(91))); end
IEVREG = round(h(92));
if (h(92) ~= -12345 & nargin == 1); disp(sprintf('IEVREG     = %d',h(92))); end
IEVTYP = round(h(93));
if (h(93) ~= -12345 & nargin == 1); disp(sprintf('IEVTYP     = %d',h(93))); end
IQUAL = round(h(94));
if (h(94) ~= -12345 & nargin == 1); disp(sprintf('IQUAL      = %d',h(94))); end
ISYNTH = round(h(95));
if (h(95) ~= -12345 & nargin == 1); disp(sprintf('ISYNTH     = %d',h(95))); end
IMAGTYP = round(h(96));
if (h(96) ~= -12345 & nargin == 1); disp(sprintf('IMAGTYP    = %d',h(96))); end
IMAGSRC = round(h(97));
if (h(97) ~= -12345 & nargin == 1); disp(sprintf('IMAGSRC    = %d',h(97))); end

%read logical header variables
%---------------------------------------------------------------------------
LEVEN = round(h(106));
if (h(106) ~= -12345 & nargin == 1); disp(sprintf('LEVEN      = %d',h(106))); end
LPSPOL = round(h(107));
if (h(107) ~= -12345 & nargin == 1); disp(sprintf('LPSPOL     = %d',h(107))); end
LOVROK = round(h(108));
if (h(108) ~= -12345 & nargin == 1); disp(sprintf('LOVROK     = %d',h(108))); end
LCALDA = round(h(109));
if (h(109) ~= -12345 & nargin == 1); disp(sprintf('LCALDA     = %d',h(109))); end

%read character header variables
%---------------------------------------------------------------------------
KSTNM = char(h(111:118));
if (str2double(KSTNM) ~= -12345 & nargin == 1); disp(sprintf('KSTNM      = %s', KSTNM)); end
KEVNM = char(h(119:134));
if (str2double(KEVNM) ~= -12345 & nargin == 1); disp(sprintf('KEVNM      = %s', KEVNM)); end
KHOLE = char(h(135:142));
if (str2double(KHOLE) ~= -12345 & nargin == 1); disp(sprintf('KHOLE      = %s', KHOLE)); end
KO = char(h(143:150));
if (str2double(KO) ~= -12345 & nargin == 1); disp(sprintf('KO         = %s', KO)); end
KA = char(h(151:158));
if (str2double(KA) ~= -12345 & nargin == 1); disp(sprintf('KA         = %s', KA)); end
KT0 = char(h(159:166));
if (str2double(KT0) ~= -12345 & nargin == 1); disp(sprintf('KT0        = %s', KT0)); end
KT1 = char(h(167:174));
if (str2double(KT1) ~= -12345 & nargin == 1); disp(sprintf('KT1        = %s', KT1)); end
KT2 = char(h(175:182));
if (str2double(KT2) ~= -12345 & nargin == 1); disp(sprintf('KT2        = %s', KT2)); end
KT3 = char(h(183:190));
if (str2double(KT3) ~= -12345 & nargin == 1); disp(sprintf('KT3        = %s', KT3)); end
KT4 = char(h(191:198));
if (str2double(KT4) ~= -12345 & nargin == 1); disp(sprintf('KT4        = %s', KT4)); end
KT5 = char(h(199:206));
if (str2double(KT5) ~= -12345 & nargin == 1); disp(sprintf('KT5        = %s', KT5)); end
KT6 = char(h(207:214));
if (str2double(KT6) ~= -12345 & nargin == 1); disp(sprintf('KT6        = %s', KT6)); end
KT7 = char(h(215:222));
if (str2double(KT7) ~= -12345 & nargin == 1); disp(sprintf('KT7        = %s', KT7)); end
KT8 = char(h(223:230));
if (str2double(KT8) ~= -12345 & nargin == 1); disp(sprintf('KT8        = %s', KT8)); end
KT9 = char(h(231:238));
if (str2double(KT9) ~= -12345 & nargin == 1); disp(sprintf('KT9        = %s', KT9)); end
KF = char(h(239:246));
if (str2double(KF) ~= -12345 & nargin == 1); disp(sprintf('KF         = %s', KF)); end
KUSER0 = char(h(247:254));
if (str2double(KUSER0) ~= -12345 & nargin == 1); disp(sprintf('KUSER0     = %s', KUSER0)); end
KUSER1 = char(h(255:262));
if (str2double(KUSER1) ~= -12345 & nargin == 1); disp(sprintf('KUSER1     = %s', KUSER1)); end
KUSER2 = char(h(263:270));
if (str2double(KUSER2) ~= -12345 & nargin == 1); disp(sprintf('KUSER2     = %s', KUSER2)); end
KCMPNM = char(h(271:278));
if (str2double(KCMPNM) ~= -12345 & nargin == 1); disp(sprintf('KCMPNM     = %s', KCMPNM)); end
KNETWK = char(h(279:286));
if (str2double(KNETWK) ~= -12345 & nargin == 1); disp(sprintf('KNETWK     = %s', KNETWK)); end
KDATRD = char(h(287:294));
if (str2double(KDATRD) ~= -12345 & nargin == 1); disp(sprintf('KDATRD     = %s', KDATRD)); end
KINST = char(h(295:302));
if (str2double(KINST) ~= -12345 & nargin == 1); disp(sprintf('KINST      = %s', KINST)); end


if nargin > 1
  for nrecs = 1:(nargin-1);
    varargout{nrecs} = eval(varargin{nrecs});
  end
end
