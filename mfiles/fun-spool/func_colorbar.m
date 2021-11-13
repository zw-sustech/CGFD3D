function [h0]=func_colorbar(n,varargin) %varargout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%              Matlab　colorbar　introduction                          %
%%%%%%%%%%%%%%%%%%%%%Non_linear colorbar%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% programmer:  w5400@pku.edu.cn 　　　　　　　　　　　　　　　　　　  %
%% @@This programe to use my own colorbar table to pic and can pic  the%
%% the picture ungradually or nonlinear@@@                             %
%%---------------------------------------------------------------------%
%%a is a 64*3 array include 64 color                                   %
%%a_n is my own non gradual color table                                %
%%colorbar has many property include 'pos','XTick,'Label','XTickLabel'.%
%%n is the band of colorbar ; pos is [Xmin,Ymin,dx,dy] of colorbar     %
%%Matlab has it's own command to do this ungradually colorbar.         %
%%eg1:colormap(jet(11));colorbar                                       %
%%eg2:contourf(peaks(30),[-8:2:6])                                     %
%%   contourcmap(1,'jet','YLabelString','Value','FontSize',6) ;colobar %
%%Matlab has it's own command to do this nonlinear colorbar.           %
%%eg: contourfm(x,y,z,[-120,-90,-30,0,40,50,100])                      %
%%---------------------------------------------------------------------%
%function(v,'colbar','nonlinear','levels',11,'cmin',0','cmax',22,'dc',2)
%function(v,'colbar','nonlinear',[10,15,26,40])                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----there is a special point in matlab for color appliying. it is :---%
%%----------------------color&data introdution-------------------------%
%   Graphics objects that use pseudocolor  -- SURFACE and PATCH objects%
%   which are created by the functions MESH, SURF, and PCOLOR -- map   %
%   a color matrix, C, whose values are in the range [Cmin, Cmax],     %
%   to an array of indices, k, in the range [1, m].                    %
%   The values of Cmin and Cmax are either min(min(C)) and max(max(C)),%
%   or are specified by CAXIS.  The mapping is linear, with Cmin       %
%   mapping to index 1 and Cmax mapping to index m.  The indices are   %
%   then used with the colormap to determine the color associated      %
%   with each matrix element.  See CAXIS for details.                  %
%%---------------------------------------------------------------------%
%caxis([min,max])is correspond to the [min,max]of colortable and meanly%
%_correspond to the data value form the min to the max.==the colortable%  
%is linear to the data range.we can't define the colorbar non_linear   %
%except the min and max.                                               %                                                 
% and if you define a caxis ,then the value>caxis_max will be shown on %
% caxis_max.the value <caxis_min will be shown on caxis_min.           %
%the function need cmin(min of caxis) cmax dc(the step length of caxis%%
%----------------------------------------------------------------------%
%if  you want to change a pic color table                              %
%only thing you must do is ,graph the pic out first in default colormap%
%then ,change the colormap . now ,the pic will change auto.            %
%this change has no connection with having a colorbar or not.          %
%eg: surf(peaks(20)), colormap cool ,( colorbar)                       %
%varargin:how to use it? eg:func_colorbarNO2(v,11,'cmin',11,'cmax',32')%
%the band of a_n must = band of colorbar;colormap(a_n(11-n+1:1:11,:))  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 if 1
a_n=[ 0 0 255 
      0 153 255
      0 255 255
      0 255 0
      128 255 0
      255 255 255
      255 200 0
      255 153 0
      255 0 0
      255 0 255
      153 0 153 ]/255

   a_n=[ 255 255 255
         200 255 255
         0 255 255
         0 153 255
         128 255 0
         0 255 0
         255 255 0
         255 153 0
         255 0 0
         255 0 255
         153 0 153 ]/255;
 end

 a_n = [          0         0    0.5625
         0         0    0.9375
         0    0.3125    1.0000
         0    0.6875    1.0000
    0.0625    1.0000    0.9375
    0.4375    1.0000    0.5625
    0.8125    1.0000    0.1875
    1.0000    0.8125         0
    1.0000    0.4375         0
    1.0000    0.0625         0
    0.6875         0         0];
if 0
 a_n = [ 0         0    0.5625
         0         0    0.9375
         0    0.3125    1.0000
         0    0.6875    1.0000
    0.0625    1.0000    0.9375
    0.4375    1.0000    0.5625
0         0    0.5625
0.0625    1.0000    0.9375
0.6875         0         0
    0.8125    1.0000    0.1875
    1.0000    0.8125         0
];     
 end 
 if 1
    a_n=[218 165 32
        255	255	255
        218 112 214
        100 149 237  
%124 252 0 
50 205 50
0 255 0
0 153 0
     
255 153 0
210 105 30
 255 69 0
 255 0 0
 135 0 0
 
%199 21 133
% 255 190 0 
]./255;

end

 a_n = [          0         0    0.5625
                1 1 1
                0.8549    0.6875    0.8392
                0.8549    0.4392    0.8392
                0.06    0.3    1.0000
                0    0.6875    1.0000
                %0.0625    1.0000    0.9375
                0.4375    1.0000    0.5625
                0.8125    1.0000    0.1875
                1.0000    0.8125         0
                1.0000    0.4375         0
                1.0000    0.0625         0
                0.6875         0         0]; 

%%%%%%%%%%%%%%%%%%global&EastAsia&china&northchina%%%%%%%%%%%%%%%%%%%%%%%
cmin=0; cmax=22; dc=2; vv = 1%default
%cmin=11; cmax=33; dc=2;
%cmin=6; cmax=24; dc=3
%cmin=11,cmax=32,dc=3
%if nargin < 11, d=0.03; end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d=0.03; %nargin is the  number of input arguments 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(varargin)%length=max(size(X))
    if ( strcmpi(varargin{i},'cmin')) 
        cmin=varargin{i+1};
    end
    if ( strcmpi(varargin{i},'cmax')) 
        cmax=varargin{i+1};
    end
    if ( strcmpi(varargin{i},'dc')) 
        dc=varargin{i+1};
    end
    if ( strcmpi(varargin{i},'h')) 
        h=varargin{i+1};
    end
    %if ( strcmpi(varargin{i},'colbar')) 
        %colbar=varargin{i+1};
    if ( strcmpi(varargin{i},'colbar')) 
        colbar=1;
    end
    if ( strcmpi(varargin{i},'nonlinear')) 
        nonlinear=1;
        vv=varargin{i+1};
    end
    if ( strcmpi(varargin{i},'gradu')) 
        gradu=1;
    end
    %if length(varargin)>0 & isnumeric(varargin{1}) & length(varargin{1})>1;
    %     	vv = varargin{1}; % contour levels
    %else
    %      	vv= [];
    %end
end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch exist('nonlinear')
    case 1
        n=length(vv);cmin=1;cmax=n;caxis([cmin,cmax])
        Labels=[1:1:n];
        Labels1=vv;  
        if n>7
                colormap(a_n(11-n+1+1:1:11,:));
        elseif n<=7
                colormap(a_n(1:2:n*2-2,:));
        end     
    case 0
        if gradu==0
        dc=(cmax-cmin)/n;
        caxis([cmin,cmax])
        Labels=[cmin:dc:cmax]
        Labels1=Labels
        %Labels1(size(Labels,2))=round(max(max(v)));
        if n>6
            colormap(a_n(11-n+1:1:11,:));
        elseif n<=6
            colormap(a_n(1:2:n*2-1,:));
        end   
        end
end
if exist('colbar')
    switch h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%colorbar's dirction(h or not)%
             case 0
                 h0=colorbar('vert');%tp(4)=dy;%dy is the hight of colorbar
                 tp=get(h0,'pos')
                 tp(3)=d;%dy is the hight of colorbar
                 set(h0,'pos',tp);%[Xmin,Ymin,dX,dY]
                 if exist('nonlinear')
                 set(h0,'YTick',Labels)
                 set(h0,'YTickLabel',Labels1)
                 end
             case 1
                 h0=colorbar('h');
                 tp=get(h0,'pos')
                 tp(4)=d;             
                 set(h0,'pos',tp);
                 if exist('nonlinear')
                     set(h0,'XTick',Labels)
                     set(h0,'XTickLabel',Labels1)
                 end
         end
 end
if nargout==1,varargout{1}=h0;end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
