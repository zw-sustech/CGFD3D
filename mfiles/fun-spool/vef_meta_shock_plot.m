function [pid,tid]=vef_meta_shock_plot(sx,sy,sz,snm,color_spec)

nshock=length(sx);
for n=1:nshock
    hold on
    pid(n)=plot3(sx(n),sy(n),sz(n), ...
      'p','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10);
      %'p','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',18);
    %text(sx(n)+5,sy(n),sz(n),snm(n,:),'color','w');
end
tid=[];

