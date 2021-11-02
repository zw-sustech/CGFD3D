%FINDTAG  locate objects with specific tag
%
% FINDTAG(STR)
% Locate objects with specific tag STR
%
% FINDTAG(OBJECT_HANDLES,STR) 
% Restricts the search to objects listed in objhandles and their descendants.
%
% FINDTAG(...,'-depth',d)
% The depth argument d controls how many levels under the handles in objhandles 
% are traversed. Specify d as inf to get the default behavior of all levels.
% Specify d as 0 to restrict to the objects listed in OBJECT_HANDLES.
%
% Rem: In order for this to work properly, the object's tag field must be a string 
% containing keywords (or tags) separated by commas.
%
%
% See also
%  addtag, hastag, rmtag
%
% G.Rilling 12/2006
% gabriel.rilling@ens-lyon.fr


function objs = findtag(varargin)

if any(ishandle(varargin{1}))
  lobj = varargin{1};
  tag = regexptranslate('escape',varargin{2});
  objs = findobj(lobj,'-regexp','Tag',['(\W|^)',tag,'(\W|$)'],varargin{3:end});
else
  tag = regexptranslate('escape',varargin{1});
  objs = findobj('-regexp','Tag',['(\W|^)',tag,'(\W|$)'],varargin{2:end});
end
