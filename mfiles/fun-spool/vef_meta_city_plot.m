function [pid,tid]=vef_meta_fault_plot(cx,cy,cz,cnm,color_spec)

ncity=length(cx);

for n=1:ncity
    if cnm{n}(1)=='#', continue; end
    hold on
    pid(n)=plot3(cx(n),cy(n),cz(n), ...
      'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',4);
    %plot3(cx(n),cy(n),cz(n), ...
    %  'square','MarkerEdgeColor','w','MarkerFaceColor','w','MarkerSize',8);
    tid(n)=text(cx(n)+5,cy(n),cz(n),cnm{n},'color','w','fontsize',9);
end
