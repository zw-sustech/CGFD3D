%EMD_VISU  visualization of EMD and partial reconstructions (fine to coarse & coarse to fine)
%
% inputs :   - x : analyzed signal, if x is complex cemd_visu is called
%            - t : time instants
%            - imf : output of emd.m
%            - i (optional) : figure number for display
%
% outputs :  - f2c : fine to coarse reconstruction
%            - c2f : coarse to fine reconstruction
%
%
% See also
%  cemd_visu (visualization of bivariate EMD)
%  emd (EMD and bivariate EMD)
%  emdc, emdc_fix (fast implementations of EMD)
%  cemdc, cemdc_fix, cemdc2, cemdc2_fix (fast implementations of bivariate EMD)
%
% P. Flandrin, Mar. 13, 2003
% G. Rilling, last modification 3.2006
% gabriel.rilling@ens-lyon.fr

function varargout = emd_visu(x,t,imf,i);

if sum(size(t)>1)>1
  imf = t;
  t = 1:length(x);
  if(nargin==3)
    fignum = i;
  end
else
  if(nargin==4)
    fignum = i;
  end
end

if ~isreal(x)||~isreal(imf)
  if exist('fignum','var')
    cemd_visu(x,t,imf,fignum);
  else
    cemd_visu(x,t,imf);
  end
  return;
end

if length(x)~= length(t)
  error('X and T must have the same length');
end
if size(imf,2) ~= length(x)
  error('the number of columns in IMF must equal the length of X')
end


if exist('fignum','var')
  figure(fignum);
else
  figure;
end

mx = min(x);
Mx = max(x);

s = size(imf);
k = s(1);

M = max(max(abs(imf(1:k-1,:))));

subplot(k+1,1,1)
plot(t,x)
axis([t(1) t(s(2)) mx Mx])
set(gca,'YTick',[])
set(gca,'XTick',[])
ylabel(['signal'])

for j = 1:k-1
  subplot(k+1,1,j+1)
  plot(t,imf(j,:))
  axis([t(1) t(s(2)) -M M])
  set(gca,'YTick',[])
  set(gca,'XTick',[])
  ylabel(['imf',int2str(j)])
end
subplot(k+1,1,1)
title('Empirical Mode Decomposition')

subplot(k+1,1,k+1)
plot(t,imf(k,:),'r')
axis('tight')
set(gca,'YTick',[])
set(gca,'XTick',[])
ylabel('res.')

f2c = [];
f2c(1,:) = imf(1,:);

c2f = [];
c2f(1,:) = imf(k,:);

for j = 2:k
  f2c(j,:) = f2c(j-1,:) + imf(j,:);
  c2f(j,:) = c2f(j-1,:) + imf(k+1-j,:);
end

if(nargin==4)
  figure(i+1)
else
  figure
end

mx = min(x);
Mx = max(x);

s = size(f2c);
k = s(1);

M = max(max(abs(f2c(1:k-1,:))));

subplot(k+1,1,1)
plot(t,x)
axis([t(1) t(s(2)) mx Mx])
set(gca,'YTick',[])
set(gca,'XTick',[])
ylabel(['signal'])

for j = 1:k-1
  subplot(k+1,1,j+1)
  plot(t,f2c(j,:))
  axis([t(1) t(s(2)) -M M])
  set(gca,'YTick',[])
  set(gca,'XTick',[])
  ylabel(['f2c',int2str(j)])
end
subplot(k+1,1,1)
title('f2c')

subplot(k+1,1,k+1)
plot(t,f2c(k,:),'r')
mr = min(f2c(k,:));
Mr = max(f2c(k,:));
axis([t(1) t(s(2)) mr Mr])
set(gca,'YTick',[])
set(gca,'XTick',[])
ylabel('signal')

if(nargin==4)
  figure(i+2)
else
  figure
end

mx = min(x);
Mx = max(x);

s = size(c2f);
k = s(1);

M = max(max(abs(c2f(1:k-1,:)-mean(x))));

subplot(k+1,1,1)
plot(t,x)
axis([t(1) t(s(2)) mx Mx])
set(gca,'YTick',[])
set(gca,'XTick',[])
ylabel(['signal'])

for j = 1:k-1
  subplot(k+1,1,j+1)
  plot(t,c2f(j,:))
  axis([t(1) t(s(2)) -M M])
  set(gca,'YTick',[])
  set(gca,'XTick',[])
  ylabel(['c2f',int2str(j)])
end
subplot(k+1,1,1)
title('c2f')

subplot(k+1,1,k+1)
plot(t,c2f(k,:),'r')
mr = min(c2f(k,:));
Mr = max(c2f(k,:));
axis([t(1) t(s(2)) mr Mr])
set(gca,'YTick',[])
set(gca,'XTick',[])
ylabel('signal')

if nargout > 0
  varargout = {f2c,c2f};
end
