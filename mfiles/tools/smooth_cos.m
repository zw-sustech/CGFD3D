clear all;

dz1 = 5
dz2 = 6
N   = 100
L   = N * (dz1+dz2)/2;

deri = (1 + sin( -pi/2 + pi/N * [0:N]) )/2;
A = L / sum(deri);

z(1) = 0;
for n = 2 : N+1
    z(n) = z(n-1) + A*deri(n);
end
