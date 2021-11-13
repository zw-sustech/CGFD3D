%ADDTAG  add a tag to an object
%
% ADDTAG(OBJ_HANDLE,STR)
% 
% Adds the tag STR to the object referrenced by OBJ_HANDLE
% When OBJ_HANDLE is an array of handles, STR is added to all corresponding
% objects.
% 
% Rem: In order for this to work properly, the object's tag field must be a string 
% containing keywords (or tags) separated by commas.
%
%
% See also
%  hastag, findtag, rmtag
%
% G.Rilling 12/2006
% gabriel.rilling@ens-lyon.fr

function addtag(obj,str)
inds = ~hastag(obj,str);
tags = get(obj(inds),'Tag');
if ~iscell(tags)
  tags = {tags};
end
inds = find(inds);
new_tags = cellfun(@(x)[x,',',str],tags,'UniformOutput',false);
new_tags = cellfun(@(x)regexprep(x,'^,',''),new_tags,'UniformOutput',false);
arrayfun(@(ind)set(obj(inds(ind)),'Tag',new_tags{ind}),1:length(inds));
