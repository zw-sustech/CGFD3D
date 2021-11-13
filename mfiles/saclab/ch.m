%CH    change SAC header
%
%    Change SAC header variables for 
%    SAC files read in to matlab with rsac.m
%    
%    Examples:
%
%    To change the SAC variable DELTA from station KATH to
%    the matlab variable dt:
%
%    KATH=ch(KATH,'DELTA',dt);
%
%    To change the SAC variables STLA and STLO from station KATH
%    to the matlab variables lat and lon:
%
%    KATH=ch(KATH,'STLA',lat,'STLO',lon)
%
%    To change the SAC variable KT0 from station SKS to sSKS
%    for station KATH:
%
%    KATH=ch(KATH,'KT0','sSKS'); 
%
%    Note:  this program can only handle one seismogram
%    file at a time, but can change multiple header values.
%
%    by Michael Thorne (4/2004)  mthorne@asu.edu
%
%    See also:  RSAC, LH, BSAC, WSAC 

function [output]=ch(file,varargin);

% first test to see if the file is indeed a sacfile
%---------------------------------------------------------------------------
if (file(303,3)~=77 & file(304,3)~=73 & file(305,3)~=75 & file(306,3)~=69)
  error('Specified Variable is not in SAC format ...')
elseif nargin <= 1
  error('Not enough input arguments ...')
end

output = file;

for i=1:2:(nargin-1);
input = double(varargin{i});

if size(varargin{i} < 8)
  insize = length(varargin{i});
  for kk=(insize+1):8;
   input(kk) = 32; 
  end
end


if input == double('DELTA   '); output(1,3) = varargin{i+1}; end
if input == double('DEPMIN  '); output(2,3) = varargin{i+1}; end
if input == double('DEPMAX  '); output(3,3) = varargin{i+1}; end
if input == double('SCALE   '); output(4,3) = varargin{i+1}; end
if input == double('ODELTA  '); output(5,3) = varargin{i+1}; end
if input == double('B       '); output(6,3) = varargin{i+1}; end
if input == double('E       '); output(7,3) = varargin{i+1}; end
if input == double('O       '); output(8,3) = varargin{i+1}; end
if input == double('A       '); output(9,3) = varargin{i+1}; end
if input == double('T0      '); output(11,3) = varargin{i+1}; end
if input == double('T1      '); output(12,3) = varargin{i+1}; end
if input == double('T2      '); output(13,3) = varargin{i+1}; end
if input == double('T3      '); output(14,3) = varargin{i+1}; end
if input == double('T4      '); output(15,3) = varargin{i+1}; end
if input == double('T5      '); output(16,3) = varargin{i+1}; end
if input == double('T6      '); output(17,3) = varargin{i+1}; end
if input == double('T7      '); output(18,3) = varargin{i+1}; end
if input == double('T8      '); output(19,3) = varargin{i+1}; end
if input == double('T9      '); output(20,3) = varargin{i+1}; end
if input == double('F       '); output(21,3) = varargin{i+1}; end
if input == double('RESP0   '); output(22,3) = varargin{i+1}; end
if input == double('RESP1   '); output(23,3) = varargin{i+1}; end
if input == double('RESP2   '); output(24,3) = varargin{i+1}; end
if input == double('RESP3   '); output(25,3) = varargin{i+1}; end
if input == double('RESP4   '); output(26,3) = varargin{i+1}; end
if input == double('RESP5   '); output(27,3) = varargin{i+1}; end
if input == double('RESP6   '); output(28,3) = varargin{i+1}; end
if input == double('RESP7   '); output(29,3) = varargin{i+1}; end
if input == double('RESP8   '); output(30,3) = varargin{i+1}; end
if input == double('RESP9   '); output(31,3) = varargin{i+1}; end
if input == double('STLA    '); output(32,3) = varargin{i+1}; end
if input == double('STLO    '); output(33,3) = varargin{i+1}; end
if input == double('STEL    '); output(34,3) = varargin{i+1}; end
if input == double('STDP    '); output(35,3) = varargin{i+1}; end
if input == double('EVLA    '); output(36,3) = varargin{i+1}; end
if input == double('EVLO    '); output(37,3) = varargin{i+1}; end
if input == double('EVEL    '); output(38,3) = varargin{i+1}; end
if input == double('EVDP    '); output(39,3) = varargin{i+1}; end
if input == double('MAG     '); output(40,3) = varargin{i+1}; end
if input == double('USER0   '); output(41,3) = varargin{i+1}; end
if input == double('USER1   '); output(42,3) = varargin{i+1}; end
if input == double('USER2   '); output(43,3) = varargin{i+1}; end
if input == double('USER3   '); output(44,3) = varargin{i+1}; end
if input == double('USER4   '); output(45,3) = varargin{i+1}; end
if input == double('USER5   '); output(46,3) = varargin{i+1}; end
if input == double('USER6   '); output(47,3) = varargin{i+1}; end
if input == double('USER7   '); output(48,3) = varargin{i+1}; end
if input == double('USER8   '); output(49,3) = varargin{i+1}; end
if input == double('USER9   '); output(50,3) = varargin{i+1}; end
if input == double('DIST    '); output(51,3) = varargin{i+1}; end
if input == double('AZ      '); output(52,3) = varargin{i+1}; end
if input == double('BAZ     '); output(53,3) = varargin{i+1}; end
if input == double('GCARC   '); output(54,3) = varargin{i+1}; end
if input == double('DEPMEN  '); output(57,3) = varargin{i+1}; end
if input == double('CMPAZ   '); output(58,3) = varargin{i+1}; end
if input == double('CMPINC  '); output(59,3) = varargin{i+1}; end
if input == double('XMINIMUM'); output(60,3) = varargin{i+1}; end
if input == double('XMAXIMUM'); output(61,3) = varargin{i+1}; end
if input == double('YMINIMUM'); output(62,3) = varargin{i+1}; end
if input == double('YMAXIMUM'); output(63,3) = varargin{i+1}; end

if input == double('NZYEAR  '); output(71,3) = varargin{i+1}; end
if input == double('NZJDAY  '); output(72,3) = varargin{i+1}; end
if input == double('NZHOUR  '); output(73,3) = varargin{i+1}; end
if input == double('NZMIN   '); output(74,3) = varargin{i+1}; end
if input == double('NZSEC   '); output(75,3) = varargin{i+1}; end
if input == double('NZMSEC  '); output(76,3) = varargin{i+1}; end
if input == double('NVHDR   '); output(77,3) = varargin{i+1}; end
if input == double('NORID   '); output(78,3) = varargin{i+1}; end
if input == double('NEVID   '); output(79,3) = varargin{i+1}; end
if input == double('NPTS    '); output(80,3) = varargin{i+1}; end
if input == double('NWFID   '); output(82,3) = varargin{i+1}; end
if input == double('NXSIZE  '); output(83,3) = varargin{i+1}; end
if input == double('NYSIZE  '); output(84,3) = varargin{i+1}; end
if input == double('IFTYPE  '); output(86,3) = varargin{i+1}; end
if input == double('IDEP    '); output(87,3) = varargin{i+1}; end
if input == double('IZTYPE  '); output(88,3) = varargin{i+1}; end
if input == double('IINST   '); output(90,3) = varargin{i+1}; end
if input == double('ISTREG  '); output(91,3) = varargin{i+1}; end
if input == double('IEVREG  '); output(92,3) = varargin{i+1}; end
if input == double('IEVTYP  '); output(93,3) = varargin{i+1}; end
if input == double('IQUAL   '); output(94,3) = varargin{i+1}; end
if input == double('ISYNTH  '); output(95,3) = varargin{i+1}; end
if input == double('IMAGTYP '); output(96,3) = varargin{i+1}; end
if input == double('IMAGSRC '); output(97,3) = varargin{i+1}; end

if input == double('LEVEN   '); output(106,3) = varargin{i+1}; end
if input == double('LPSPOL  '); output(107,3) = varargin{i+1}; end
if input == double('LOVROK  '); output(108,3) = varargin{i+1}; end
if input == double('LCALDA  '); output(109,3) = varargin{i+1}; end

if input == double('KSTNM   ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(111:118,3) = newname(1:8)';
end
if input == double('KEVNM   ');
 newname = double(varargin{i+1});
  if size(newname) < 16; newname((length(varargin{i+1})+1):16) = 32; end
 output(119:134,3) = newname(1:16)';
end
if input == double('KHOLE   ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(135:142,3) = newname(1:8)';
end
if input == double('KO      ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(143:150,3) = newname(1:8)';
end
if input == double('KA      ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(151:158,3) = newname(1:8)';
end
if input == double('KT0     ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(159:166,3) = newname(1:8)';
end
if input == double('KT1     ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(167:174,3) = newname(1:8)';
end
if input == double('KT2     ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(175:182,3) = newname(1:8)';
end
if input == double('KT3     ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(183:190,3) = newname(1:8)';
end
if input == double('KT4     ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(191:198,3) = newname(1:8)';
end
if input == double('KT5     ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(199:206,3) = newname(1:8)';
end
if input == double('KT6     ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(207:214,3) = newname(1:8)';
end
if input == double('KT7     ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(215:222,3) = newname(1:8)';
end
if input == double('KT8     ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(223:230,3) = newname(1:8)';
end
if input == double('KT9     ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(231:238,3) = newname(1:8)';
end
if input == double('KF      ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(239:246,3) = newname(1:8)';
end
if input == double('KUSER0  ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(247:254,3) = newname(1:8)';
end
if input == double('KUSER1  ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(255:262,3) = newname(1:8)';
end
if input == double('KUSER2  ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(263:270,3) = newname(1:8)';
end
if input == double('KCMPNM  ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(271:278,3) = newname(1:8)';
end
if input == double('KNETWK  ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(279:286,3) = newname(1:8)';
end
if input == double('KDATRD  ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(287:294,3) = newname(1:8)';
end
if input == double('KINST   ');
 newname = double(varargin{i+1});
  if size(newname) < 8; newname((length(varargin{i+1})+1):8) = 32; end
 output(295:302,3) = newname(1:8)';
end

end
