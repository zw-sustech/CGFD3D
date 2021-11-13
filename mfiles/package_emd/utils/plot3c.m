%PLOT3C  plots a complex signal in 3D
%
% use: same as plot
%
%
% See also
%  plotc (another type of visualization of complex signals)
%
% G. Rilling, last modification 12.2006
% gabriel.rilling@ens-lyon.fr
function varargout = plot3c(varargin)

% if first arument is a handle to some axes object
if nargin >=2 && isscalar(varargin{1}) && ishandle(varargin{1})
  ax = varargin{1};
  varargin = varargin(2:end);
else
  ax = gca;
end

char_inds = find(cellfun(@ischar,varargin));
if ~isempty(char_inds)
  numeric_args = char_inds(1)-1;
else
  numeric_args = length(varargin);
end

switch numeric_args
  case 1
    y = varargin{1};
    x = 1:length(y);
  case 2
    y = varargin{2};
    x = varargin{1};
end
h=plot3(ax,x,real(y),imag(y),varargin{numeric_args+1:end});
if nargout
  varargout = {h};
end