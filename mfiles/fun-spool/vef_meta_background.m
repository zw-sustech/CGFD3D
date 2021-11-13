function [aid,sid]=vef_meta_backgroud(fid,x,y,z,scl_daspect,color_spec,bgcolor)

C=zeros(size(z));
sid=surf(x,y,z,C);
view(0,90)
if exist('scl_daspect')
   daspect(scl_daspect);
end
axis tight
shading flat

%c1=name2rgb('red'); c2=name2rgb('red');
%cMap=makeColorMap(c1,c2,5);
%cMap=[1 0 0; 1 0  0; 1 0 0];
cMap=[color_spec; color_spec; color_spec];
colormap(cMap);

if exist('bgcolor','var')
   set(gcf,'color',bgcolor);
end

aid=gca;
