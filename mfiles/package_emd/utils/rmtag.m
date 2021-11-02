%RMTAG  removes tag from object
%
% RMTAG(OBJ_HANDLE,STR)
% 
% Removes the tag STR from the object referrenced by OBJ_HANDLE
% When OBJ_HANDLE is an array of handles, the tag is removed from all
% corresponding objects.
% 
% Rem: In order for this to work properly, the object's tag field must be a string 
% containing keywords (or tags) separated by commas.
%
%
% See also
%  addtag, hastag, findtag
%
% G.Rilling 12/2006
% gabriel.rilling@ens-lyon.fr

function rmtag(obj,str)
if any(~hastag(obj,str))
    warning('rmtag:warning','no such tag in object')
end

arrayfun(@rmtag1,1:length(obj));

  function rmtag1(ind)
    tag = get(obj(ind),'Tag');
    tag = regexprep(tag,['(\W|^)',str,'(\W|$)'],'$1$2');
    tag = regexprep(tag,',,',',');
    tag = regexprep(tag,'^,|,$','');
    set(obj(ind),'Tag',tag);
  end

end

