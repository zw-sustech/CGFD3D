%CEMD_VISU  visualization of bivariate/complex EMD
%
% inputs :   - x : complex analyzed signal
%            - t : time instants
%            - imf : set of complex IMFs
%            - i (optional) : figure number for display
%
%
% The slider on the right hand side of the figure controls the direction on
% which the complex signal is projected:
% phi=0: the projection corresponds to the real part of the signal,
% phi=pi/2: the projection corresponds to the imaginary part of the signal.
%
% examples:
%
%>>cemd_visu(x,imf)
%>>cemd_visu(x,t,imf)
%>>cemd_visu(x,t,imf,1)
%
%
% See also
%  emd (EMD and bivariate EMD)
%  cemdc, cemdc_fix, cemdc2, cemdc2_fix (fast implementations of bivariate EMD)
%
% G. Rilling, last modification 3.2007
% gabriel.rilling@ens-lyon.fr

function cemd_visu(x,t,imf,numfig)

if sum(size(t)>1)>1
  imf = t;
  t = 1:length(x);

  if(nargin==3)
    figure(numfig)
  else
    figure
  end
else
  if(nargin==4)
    figure(numfig)
  else
    figure
  end

end

if length(x)~= length(t)
  error('X and T must have the same length');
end
if size(imf,2) ~= length(x)
  error('the number of columns in IMF must equal the length of X')
end


set(gcf,'name','cemd_visu: EMD','Toolbar','figure')

phi = 0;

[k,n] = size(imf);

mx = mean(x);
Mx = max(abs(x-mx));
mr = mean(imf(end,:));
Mr = max(abs(imf(end,:)-mr));
M = max(max(abs(imf(1:k-1,:))));
A_x = [t(1),t(n),real(mx)-Mx,real(mx)+Mx];
A_imf = [t(1),t(n),-M,M];
A_res = [t(1),t(n),real(mr)-Mr,real(mr)+Mr];

ax_x = subplot(k+1,1,1);
set(gca,'NextPlot','replacechildren','YTick',[],'XTickLabel',{},'Box','on')
axis(A_x);
ylabel(['signal'])
T=title('Empirical Mode Decomposition, \phi/2\pi=0');

for j=1:k-1
  ax_imf(j) = subplot(k+1,1,j+1);
  set(gca,'NextPlot','replacechildren','YTick',[],'XTickLabel',{},'Box','on')
  ylabel(['imf',int2str(j)])
  axis(A_imf);
end
ax_res = subplot(k+1,1,k+1);
set(gca,'NextPlot','replacechildren','Box','on')
axis(A_res);
ylabel('res.')

plot_emd(phi);

P(1,:) = get(ax_x,'Position');
P(2,:) = get(ax_res,'Position');
L = 1-(P(1,1)+P(1,3));
H = P(1,2)+P(1,4)-P(2,2);
slider_handle = uicontrol('Style','Slider','Min',0,'Max',1,'Value',0,'SliderStep',[.01,.05],'Units','normalized','Position',[1-2*L/3,P(2,2),L/3,H],'Callback',@slider_callback);

  function plot_emd(phi)

    plot(ax_x,t,real(exp(-i*phi)*x));
    set(ax_x,'Ylim',[real(exp(-i*phi)*mx)-Mx,real(exp(-i*phi)*mx)+Mx]);

    for j = 1:k-1
      plot(ax_imf(j),t,real(exp(-i*phi)*imf(j,:)));
    end
    plot(ax_res,t,real(exp(-i*phi)*imf(end,:)),'r');
    set(ax_res,'Ylim',[real(exp(-i*phi)*mr)-Mr,real(exp(-i*phi)*mr)+Mr]);
    set(T,'String',['Empirical Mode Decomposition, \phi/2\pi=',num2str(phi/(2*pi))]);
    drawnow;
  end

  function slider_callback(varargin)
    v = get(slider_handle,'Value');
    phi = 2*pi*v;
    plot_emd(phi);
  end

end
