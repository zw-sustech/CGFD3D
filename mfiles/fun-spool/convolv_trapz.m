function [W]=convolv_trapz(U,V,dt,nt1,nt2)
%
% convolv_trapz: Convolution using 2th order trapz.
%
% Usage: [W]=convolv_trapz(U,V,dt,nt1,nt2)
%

% Major ChangeLog:
%   2009-01-09 Wei Zhang
%     * Added help information, but uncomplete.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Date$
% $Revision$
% $LastChangedBy$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W=zeros(length(U),1);
for i=nt1:nt2
   W(i)=0.5*dt*U(1)*V(i)+0.5*dt*U(i)*V(1);
   for j=2:i-1
      W(i)=W(i)+dt*U(j)*V(i-j+1);
   end
end

%---------------------------------------------------------------------
%function [W]=cal_convolv(U,V,dt,nt1,nt2)
%nu=length(U);
%Y=conv(U,V);
%W=Y(1:nu)*dt; % U(1),V(1)=0, so equivalent

%function [W]=cal_convolv(U,V,dt,nt1,nt2)
%nu=length(U);
%W(1:nu)=0.0;
%for i=nt1:nt2
%    W(i)=sum(sum(U(1:i).*V(i:-1:1)))*dt;
%end
