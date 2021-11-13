% TRIANG.M
%
% P. Flandrin, Mar. 13, 2003
%
% generates a triangular waveform
%
% inputs :   - N : # of data samples
%            - p : period
%
% output :   - x : signal

function x = triang(N,p);

K = ceil(N/(2*p-2));

w = zeros(1,K*(2*p-2)+1);
rp = linspace(-1,1,p);
rn = linspace(1,-1,p);
r = [rp rn(2:p)];

for k = 1:K
	w(1+(k-1)*(2*p-2):1+k*(2*p-2)) = r;
end

x = w(1:N);
