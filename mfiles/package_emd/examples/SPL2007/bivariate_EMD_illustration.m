%bivariate_EMD_illustration.m 
%illustration of the bivariate EMD extension on a real-world oceanographic signal
%reproduces Fig. 3 in "Bivariate Empirical Mode Decomposition", G. Rilling,
%P. Flandrin, P. Goncalves and J. M. Lilly, IEEE Signal Processing Letters
%
%G. Rilling 3/2007 email:  gabriel.rilling@ens-lyon.fr

load('float_position_record.mat','x');

[imf,nb] = cemdc2_fix([],x,10,[],32);
n = size(imf,1);

figtitle1 = 'Float position record';
figure('name',figtitle1)
plot(x);
xlabel('Displacement East (km) --- Real part')
ylabel('Displacement North (km) --- Imaginary part')
title(figtitle1)
axis equal;
set(gca,'Ylim',[-250,300])

figtitle2 = 'Bivariate Empirical Mode Decomposition of Float signal';
figure('name',figtitle2)
subplot(n+1,1,1)
plot(real(x))
hold on
plot(imag(x),'k--')
axis tight
ylabel('signal')
title(figtitle2)
set(gca,'XTickLabel',{})
minmin = @(x)min(x(:));
maxmax = @(x)max(x(:));
m = minmin([real(imf(1:end-1,:));imag(imf(1:end-1,:))]);
M = maxmax([real(imf(1:end-1,:));imag(imf(1:end-1,:))]);
for k = 1:n
  subplot(n+1,1,k+1)
  plot(real(imf(k,:)))
  hold on
  plot(imag(imf(k,:)),'k--')
  axis([1,length(x),m,M])
  ylabel(['d_',int2str(k)])
  if k<n
    set(gca,'XTickLabel',{})
  end
end
ylabel('res.')
xlabel('Time (days)')
axis tight
