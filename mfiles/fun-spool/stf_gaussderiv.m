function [S]=fun_gaussderiv(t,a,t0);
S=zeros(size(t));
indx=find(t>=0.0 & t<=2*t0);
S(indx)=exp(-(t(indx)-t0).^2./(a^2))./(sqrt(pi)*a).*(-2.*(t(indx)-t0)./a^2);

