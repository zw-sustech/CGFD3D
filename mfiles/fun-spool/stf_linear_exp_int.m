function S=stf_linear_exp_int(T,t0)

nt=length(T);
S=zeros(nt,1);

S = 1-(1+T/t0).*exp(-T/t0);

