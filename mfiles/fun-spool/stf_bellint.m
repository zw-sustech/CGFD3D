function S=stf_bellint(T,t0)

nt=length(T);
S=zeros(nt,1);

indx=find(T>=0.0 & T<=t0);
S(indx)=T(indx)/t0 - 1/(2*pi) * sin(2*pi*T(indx)/t0);

indx=find(T>=t0);
S(indx)=1;


