function lid=vef_light(varargin)

view(-20,35);
set(gca,'box','off');
%if ~exist('lid','var'),
   lid=camlight('headlight','local');
%end
camlight(lid,70,8,'local')
%camlight(lid,-30,60)
lighting phong
set(gca,'box','on');
%lid=camlight('headlight','infinite')
%camlight(lid,40,15,'local')
%camlight(lid,30,40) 
