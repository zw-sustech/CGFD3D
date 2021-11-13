%PLOTC  plots a complex signal in 2D projection with variable projection angle.
%The angle can be modified at any time using the slider:
% slider at bottom: angle=0 the projection is the real part of the signal
% slider at 1/4: angle=pi/2 the projection is the imaginary part of the signal
% slider at top: angle=2pi same as 0
%
%Additionally, the user can choose to lock the axes or not when the angle
%is changed through the axes context menu. Note that the latter is disabled
%when tools from the figure toolbar (zoom,...) are selected.
%
% use: same as plot
%
% rem: multiple uses of plotc in the same axes are possible: it
%      produces only one slider that controls all the plots.
%
% See also
%  plot3c (another type of visualization of complex signals)
%
% G. Rilling, last modification 3.2007
% gabriel.rilling@ens-lyon.fr
function varargout = plotc(varargin)

slider_prop = 0.05;
slider_dec = 0.01;
lock_axes = 0;

% input processing
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
    if isvector(y)
      x = 1:length(y);
      y = y(:);
    else
      x = 1:size(y,1);
    end
  case 2
    y = varargin{2};
    x = varargin{1};
end
args = varargin(numeric_args+1:end);

fig = ancestor(ax,'Figure');

% context menu and slider
cmenu = get(ax,'UIContextMenu');
if isempty(cmenu)
  cmenu = uicontextmenu;
  set(ax,'UIContextMenu',cmenu);
end
if ~hastag(ax,'complexplot')
  lock_axes_menu_item = uimenu(cmenu,'Label','Lock axes','Callback',@togglelock);
  P = get(ax,'Position');
  slider = uicontrol('Style','Slider','Min',0,'Max',1,'Value',0,'SliderStep',[.01,.05],'Units','normalized','Position',[P(1)+P(3)+P(3)*slider_dec,P(2),P(3)*slider_prop,P(4)],'Callback',@slider_callback);
  addtag(slider,'complexplotslider')
else
  slider = findtag(fig,'complexplotslider');
end

% setup the complex plot
phi = get_axes_phase;
obj = plot(ax,x,real(exp(-i*phi)*y),args{:});
addtag(ax,'complexplot') % needs to be after the plot command in case the latter redefines the Tag property of the axes (default behavior)
addtag(obj,'complex');
arrayfun(@(x)setappdata(obj(x),'complex_data',y(:,x)),1:length(obj));
set(obj,'UIContextMenu',cmenu);
if nargout
  varargout = {obj};
end

set(fig,'toolbar','figure');

  function show_phase(phi)
    Xl = get(ax,'Xlim');
    Yl = get(ax,'Ylim');
    ax_ch = findtag(ax,'complex');
    for child = ax_ch'
      z = getappdata(child,'complex_data');
      z = real(exp(-i*phi)*z);
      set(child,'YData',z);
    end
    set(ax,'Xlim',Xl);
    if ~lock_axes
      set(ax,'YlimMode','auto');
    else
      set(ax,'Ylim',Yl);
    end
    drawnow;
  end

  function slider_callback(varargin)
    show_phase(get_axes_phase());
  end

  function phi = get_axes_phase()
    v = get(slider,'Value');
    phi = 2*pi*v;
  end

  function togglelock(varargin)
    lock_axes = ~lock_axes;
    if lock_axes
      set(lock_axes_menu_item,'Checked','on');
    else
      set(lock_axes_menu_item,'Checked','off');
    end
  end

end
