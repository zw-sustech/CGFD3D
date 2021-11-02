%WSAC    Write SAC binary files.
%    WSAC('sacfile') writes a SAC (seismic analysis code) binary
%    format file 
%
%    Default byte order is big-endian.  M-file can be set to default
%    little-endian byte order.
%
%    Examples:
%
%    wsac('KATH.R',kath);
%
%    wsac('SQRL.R',sqrl,'AAK.R',aak);
%
%    by Michael Thorne (4/2004)   mthorne@asu.edu
%
%    See also:  RSAC, LH, CH, BSAC 

function wsac(varargin);

for nrecs = 1:2:(nargin-1)

  outname = varargin{nrecs};
  sacfile = varargin{nrecs+1};

% first test to see if the file is indeed a sacfile
%---------------------------------------------------------------------------
if (sacfile(303,3)~=77 & sacfile(304,3)~=73 & sacfile(305,3)~=75 & ...
  sacfile(306,3)~=69)
  error('Specified Variable is not in SAC format ...')
elseif nargin <= 1
  error('Not enough input arguments ...')
end

%---------------------------------------------------------------------------
%    Default byte-order
%    endian  = 'big'  big-endian byte order (e.g., UNIX)
%            = 'lil'  little-endian byte order (e.g., LINUX)

%endian = 'big-endian';
endian = 'little-endian';

if strcmp(endian,'big-endian')
  fid = fopen(outname,'w','ieee-be'); 
elseif strcmp(endian,'little-endian')
  fid = fopen(outname,'w','ieee-le'); 
end

%read in 3-column sac-like matlab files
  h(1:306) = sacfile(1:306,3);

% write single precision real header variables:
%---------------------------------------------------------------------------
for i=1:70
  fwrite(fid,h(i),'single');
end

% write single precision integer header variables:
%---------------------------------------------------------------------------
for i=71:105
  fwrite(fid,h(i),'int32');
end

% write logical header variables
%---------------------------------------------------------------------------
for i=106:110
  fwrite(fid,h(i),'int32');
end

% write character header variables
%---------------------------------------------------------------------------
for i=111:302
  fwrite(fid,h(i),'char');
end

% write out amplitudes
%---------------------------------------------------------------------------

YARRAY = sacfile(:,2);
fwrite(fid,YARRAY,'single');

fclose(fid);

end
