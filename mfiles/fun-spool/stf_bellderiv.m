function S=stf_bellderiv(T,t0)

nt=length(T);
S=zeros(nt,1);

indx=find(T>=0.0 & T<=t0);
S(indx)=2*pi/t0^2 * sin(2*pi*T(indx)/t0);

