function [ANGRT,MDLA,MDLO,xd,yd,zd]=fun_plane4slice(xscl,yscl,zscl,STLA,STLO,EVLA,EVLO,zval)

NSAMP=1000;

MDLA=(STLA+EVLA)/2; MDLO=(STLO+EVLO)/2;
%ANGRT=atand((STLA-EVLA)/(STLO-EVLO))
ANGRT=atan2((STLA-EVLA),(STLO-EVLO))/pi*180;
if ANGRT>90, ANGRT=ANGRT-180; end
if ANGRT<-90, ANGRT=ANGRT+180; end

% plane
hspx=linspace(xscl(1),xscl(2),NSAMP);
hspy=linspace(yscl(1),yscl(2),NSAMP);
hspz=linspace(zscl(1),zscl(2),NSAMP);

[horx,hory]=meshgrid(hspx,hspy);

[matx,matz]=meshgrid(hspx,hspz);

yd=[];xd=[];zd=[];
hidpln=figure('visible','off');
% great plane
hsp=surf(matx,ones(1000)*MDLA,matz);
rotate(hsp,[0 0 1],ANGRT,[MDLO,MDLA,0]);
xd{end+1}=get(hsp,'XData'); yd{end+1}=get(hsp,'YData'); zd{end+1}=get(hsp,'ZData');
delete(hsp)

% perpendicular plane
hsp=surf(matx,ones(1000)*MDLA,matz);
rotate(hsp,[0 0 1],ANGRT+90,[MDLO,MDLA,0]);
xd{end+1}=get(hsp,'XData'); yd{end+1}=get(hsp,'YData'); zd{end+1}=get(hsp,'ZData');
delete(hsp)

% quarter plane
%hsp=surf(hory,horx,zeros(1000));
%yd{end+1}=get(hsp,'XData'); xd{end+1}=get(hsp,'YData'); zd{end+1}=get(hsp,'ZData');
%zd{end}=zd{end}+100;
%delete(hsp)

% quarter plane
hsp=surf(horx,hory,zeros(1000));
if (ANGRT>=0 & ANGRT <=90) | (ANGRT<-90)
   ANGRT4=min(abs(ANGRT),180-abs(ANGRT))
else
   ANGRT4=90-min(abs(ANGRT),180-abs(ANGRT))
end
rotate(hsp,[0 0 1],ANGRT4,[hspx(1),hspy(1),0]);
xd{end+1}=get(hsp,'XData'); yd{end+1}=get(hsp,'YData'); zd{end+1}=get(hsp,'ZData');
xd{end}=xd{end}-hspx(1)+MDLO; yd{end}=yd{end}-hspy(1)+MDLA;
zd{end}=zd{end}+zval(1);
delete(hsp)

if length(zval)>=2
% quarter plane a
ANGRT4a=ANGRT4-90;
hsp=surf(horx,hory,zeros(1000));
rotate(hsp,[0 0 1],ANGRT4a,[hspx(1),hspy(1),0]);
xd{end+1}=get(hsp,'XData'); yd{end+1}=get(hsp,'YData'); zd{end+1}=get(hsp,'ZData');
xd{end}=xd{end}-hspx(1)+MDLO; yd{end}=yd{end}-hspy(1)+MDLA;
%zd{end}=zd{end}+80;
zd{end}=zd{end}+zval(2);
delete(hsp)
end

if length(zval)>=3
% quarter plane b
ANGRT4b=ANGRT4+90;
hsp=surf(horx,hory,zeros(1000));
rotate(hsp,[0 0 1],ANGRT4b,[hspx(1),hspy(1),0]);
xd{end+1}=get(hsp,'XData'); yd{end+1}=get(hsp,'YData'); zd{end+1}=get(hsp,'ZData');
xd{end}=xd{end}-hspx(1)+MDLO; yd{end}=yd{end}-hspy(1)+MDLA;
%zd{end}=zd{end}+80;
zd{end}=zd{end}+zval(3);
delete(hsp)
end

if length(zval)>=4
% quarter plane b
ANGRT4b=ANGRT4+180;
hsp=surf(horx,hory,zeros(1000));
rotate(hsp,[0 0 1],ANGRT4b,[hspx(1),hspy(1),0]);
xd{end+1}=get(hsp,'XData'); yd{end+1}=get(hsp,'YData'); zd{end+1}=get(hsp,'ZData');
xd{end}=xd{end}-hspx(1)+MDLO; yd{end}=yd{end}-hspy(1)+MDLA;
%zd{end}=zd{end}+80;
zd{end}=zd{end}+zval(4);
delete(hsp)
end

% x 1 plane
%hsp=surf(maty,ones(1000)*MDLA,matz);
%%rotate(hsp,[0 0 1],ANGRT,[MDLO,MDLA,0]);
%yd{end+1}=get(hsp,'XData'); xd{end+1}=get(hsp,'YData'); zd{end+1}=get(hsp,'ZData');
%xd{end}=xd{end}+(MDLA-MDLA*(scl_ylim(2)-MDLO)/tan((90-ANGRT4)/180*pi));
%zd{end}=zd{end}+50;
%delete(hsp)

% y end plane
%hspxz=linspace(50,400,1000);
%[matx,matxz]=meshgrid(hspx,hspxz);
%hsp=surf(ones(1000)*scl_ylim(end),matx,matxz);
%%rotate(hsp,[0 0 1],ANGRT,[MDLO,MDLA,0]);
%yd{end+1}=get(hsp,'XData'); xd{end+1}=get(hsp,'YData'); zd{end+1}=get(hsp,'ZData');
%xd{end}=xd{end}-hspx(1)+(MDLA-(scl_ylim(2)-MDLO)*tan((90-ANGRT4)/180*pi));
%%zd{end}=zd{end}+50;
%indx=find(zd{end}>400); xd{end}(indx)=[];yd{end}(indx)=[];zd{end}(indx)=[];
%delete(hsp)

close(hidpln)

