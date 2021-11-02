%HASTAG  tests if an object has a specific tag
%
% BOOL = HASTAG(OBJ_HANDLE,STR)
% 
% Tests if the object corresponding to OBJ_HANDLE has the tag STR.
% When OBJ_HANDLE is an array of handles HASTAG returns a logical array of
% the same size.
% 
% Rem: In order for this to work properly, the object's tag field must be a string 
% containing keywords (or tags) separated by commas.
%
%
% See also
%  addtag, findtag, rmtag
%
% G.Rilling 12/2006
% gabriel.rilling@ens-lyon.fr

function bool=hastag(obj,str)
tag = get(obj,'Tag');
tmp = regexp(tag,['(\W|^)',str,'(\W|$)'],'once');
if iscell(tmp)
  bool = ~cellfun(@isempty,tmp);
else
  bool = ~isempty(tmp);
end