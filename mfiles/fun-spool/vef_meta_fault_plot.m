function [h,pid]=vef_meta_fault_plot(aid,fx,fy,fz,color_spec)

%h=newplot(aid,'NextPlot','add');
%h=newplot(aid);
hold on;
h=aid;
pid=plot3(fx,fy,fz,'k','linewidth',2,'color',color_spec);

