function seis_measure(src,event,R,S)
%fpos=get(gca,'position')
%set(gca,'position',[0.5 0.5 0.2 0.2])

fid=figure('position',[700 100 800 600], 'renderer','zbuffer','menubar','figure', ...
     'Name',['Measurement of Records ' R.name ' and Synthetic ' S.name]);
pidr=plot(R.t,R.v/R.v0,R.spec_line);
hold on
pids=plot(S.t,S.v/S.v0,S.spec_line);
%text(0.8*R.t(end),0.5*max(R.v/R.v0),R.name);
scl_ylim=get(gca,'ylim'); y0=scl_ylim(2)-(scl_ylim(2)-scl_ylim(1))/5;
text(R.t(1),y0,R.name);
