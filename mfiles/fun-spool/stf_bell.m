function S=stf_bell(T,t0)

nt=length(T);
S=zeros(nt,1);

indx=find(T>=0.0 & T<=t0);
S(indx)=(1-cos(2*pi*T(indx)/t0))/t0;

