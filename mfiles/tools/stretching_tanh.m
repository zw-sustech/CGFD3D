dz1 = 5;
dz2 = 15;
z1 = 0;
z2 = 150;
%z2 = 350;

nx  = 20;
%ds1 = 1/nx * 0.5;
%ds2 = 1/nx * 1.5;
ds1 = dz1/(z2-z1);
ds2 = dz2/(z2-z1);

A = sqrt(ds2)/sqrt(ds1);
B = 1/nx/sqrt(ds1*ds2);

x=[1:1000]*0.01;
figure
plot(x, sinh(x)./x );
%plot(x, sinh(x)./x -B);
disp('pause');
pause

delt = fzero( @(x) sinh(x)/x-B ,1e-10)
%delt = fzero( @(x) sinh(x)/x-B ,1)

for n = 1 : (nx + 1)
   u = 0.5*( 1 + tanh(delt *( (n-1)/nx - 0.5 )) / tanh(delt/2) );
   s(n) = u / ( A + (1-A)*u );
end

for n = 1 : (nx + 1)
    z(n) = z1 + (z2-z1) * s(n);
end


