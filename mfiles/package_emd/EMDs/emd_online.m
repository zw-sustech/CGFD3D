function [imf,ort,nbit] = emd_online(x,t,stop,nbpresift,tst,tst2)

%EMD_ONLINE  (On Line Empirical Mode Decomposition) computes on-line EMD
%
% IMPORTANT: EMD_ONLINE does not truly apply EMD on-line but it does AS IF.
% It is rather a demonstration that EMD can be applied on-line.
%
% stopping criterion for sifting : 
%   at each point : mean amplitude < threshold*envelope amplitude 
%                   if mean amplitude > max(envelope amplitude)/tolerance
%   &
%   at each point : mean amplitude < threshold2*enveloppe amplitude 
%                   if mean amplitude > max(envelope amplitude)/tolerance2
%
%[imf,ort,nbits] = EMD_ONLINE(x,t,stop,nbpresift,tst,tst2)
% inputs: 
%         - x : analyzed signal 
%         - t (optional) : sampling times (default : 1:length(x))
%         - stop (optional) : threshold, and threshold2 (optional)
%                             tolerance, and tolerance2 (both optional)
%                             for sifting stopping criterion 
%                             default : [0.05,0.5,20,100]
%         - nbpresift (optional) : number of sifting by blocks iterations (default 4)
%         - tst (optional) : if equals to 1 shows sifting steps
%         - tst2 (optional) : if equals to 1 shows sifting by blocks steps                   
%
% outputs:
%         - imf : intrinsic mode functions (last line = residual)
%         - ort : index of orthogonality
%         - nbits : number of iterations for each mode
%
% calls:
%         - extr (finds extrema and zero-crossings)
%         - io : computes the index of orthogonality
% 
% G. Rilling, July 2002
% gabriel.rilling@ens-lyon.fr


DEFSTOP = [0.05,0.5,20,100];% default parameters for sifting stop

NBPRESIFT = 4;%number of sifting iterations per block

if(nargin==1)
  t = 1:length(x);
  stop = DEFSTOP;
  tst = 0;
  tst2 = 0;
end

if(nargin==2)
  stop = DEFSTOP;
  tst = 0;
  tst2 = 0;
end

if (nargin==3)
  tst=0;
  tst2 = 0;
end

if (nargin==4)
  tst=0;
  tst2 = 0;
end

if (nargin==5)
  tst2 = 0;
end

if nargin > 3
    NBPRESIFT = nbpresift;
end

S = size(x);
if ((S(1) > 1) & (S(2) > 1)) | (length(S) > 2)
  error('x must have only one row or one column')
end

if S(1) > 1
  x = x';
end

S = size(t);
if ((S(1) > 1) & (S(2) > 1)) | (length(S) > 2)
  error('t must have only one row or one column')
end

if S(1) > 1
  t = t';
end

if (length(t)~=length(x))
  error('x and t must have the same length')
end

S = size(stop);
if ((S(1) > 1) & (S(2) > 1)) | (S(1) > 4) | (S(2) > 4) | (length(S) > 2)
  error('stop must have only one row or one column of max four elements')
end

if S(1) > 1
  stop = stop';
  S = size(stop);
end

if S(2) < 4
    stop(4) = DEFSTOP(4);
end

if S(2) < 3
    stop(3) = DEFSTOP(3);
end

if S(2) < 2
    stop(2) = DEFSTOP(2);
end

if S(2) == 1
  stop=[stop, DEFSTOP(2)];
end


sd = stop(1);
sd2 = stop(2);
tol = stop(3);
tol2 = stop(4);

if tst
  figure
  figures(1) = gcf;
  figure
  figures(3) = gcf;
end

if tst2
  figure
  figures(2) = gcf;
end

  
MAXITERATIONS=10000;
LARGMIN = 5;

NBSYM = 2;% maximum number of symmetrized points for interpolations

LARGTRANS = 10;
LARGTRANSPS = 5;

PAS = 20;

STEP = 5;% maximal number of iterations on a mode




LX = length(x);

% for display
sdt(LX) = 0;
sdt = sdt+sd;
sd2t(LX) = 0;
sd2t = sd2t+sd2;


% number of minima and maxima on the considered zone
lm = 0;
lM = 0;

% number of minima and maxima right of the considered zone, 
% after "stop" or "stopps"
lmr = 0;
lMr = 0;

% same, but left before "start"
lml = 0;
lMl = 0;

% total number of extrema, left and right
nem = 0;
nemr = 0;
neml = 0;

k = 1;
nbit = 0;

% number of modes, and number of modes on which block siftings are completed
nbmodes = 1;
nbmodes_psdone = 0;

start = 1;

% end of the constant part of the window
stop = min(PAS+1,LX);

% start and end of the considered zone
stopr = 1;
startl = 1;

% end of available data on the considered zone
fin = 1;

% start of the considered zone for block sifting
limpsl(1,1:NBPRESIFT) = 1;

% start and end of segment to which sifting is applied
startps = 1;
stopps(1,NBPRESIFT) = 0;
stopps = stopps + 1;

% end of available data for block sifting
finps(1,1:NBPRESIFT) = 1;
finps(1,1) = 10*PAS;

% tests if all data are available for an iteration of block sifting
lafinps(1,NBPRESIFT) = 0;

% allows to interrupt a mode extraction to process to the next mode
% interrupts also if not enough data available
suspps(1,NBPRESIFT) = 0;

% tests for the termination of one iteration of block sifting 
stoptestps(1,NBPRESIFT) = 0;

indmin = [];
indmax = [];

% tests if all data are available for on-line sifting
lafin = 0;

% tests if a mode is entirely extracted
stoptest = 0;

% allows to interrupt a mode extraction to process to the next mode
% interrupts also if not enough data available
susp = 0;

% tells if the considered zone has to be moved forward for having enough extrema
needextr = 1;

% idem for block sifting
needextrps(1,1:NBPRESIFT) = 1;

% tells how many iterations of block sifting have been initiated
nbstartedpresift = 1;

% modes concerned by block sifting
mps = x;

% mode concerned by on-line sifting
m(LX) = 0;

trig = 0;

if tst | tst2
  disp('appuyer sur une touche pour commencer')
  pause
end
  
while sum(stoptest) < nbmodes % global loop

  for k = 1:nbmodes

	  nsteps = 0;
    waittest = 0;

    if k == 1 & trig
        suspps(1,1) = 0;
        trig = 0;
    end
    
    
    while  sum(stoptestps(k,:)) < NBPRESIFT & ~waittest & sum(suspps(k,:)) < nbstartedpresift(k) % boucle de presifting
      for i = 1:nbstartedpresift(k)

          
	if needextrps(k,i) == 1
	  [indmintmp,indmaxtmp] = extr(mps(k,max([(limpsl(k,i)-1),1]):finps(k,i),i));
	  nb = sum(indmintmp > stopps(k,i))+sum(indmaxtmp > stopps(k,i));
	  stoprps(k,i) = finps(k,i);
      
	  if nb < 8*LARGTRANSPS & finps(k,i) < LX
          suspps(k,i) = 1;
          if k == 1 & i == 1
              finps(1,1) = min(LX,finps(1,1) + 10*PAS);
              [indmintmp,indmaxtmp] = extr(mps(k,max([(limpsl(k,i)-1),1]):finps(k,i),i));
              nb = sum(indmintmp > stopps(k,i))+sum(indmaxtmp > stopps(k,i));
              stoprps(k,i) = finps(k,i);
              trig = 1;
          end
      else

	    lmt = length(indmintmp);
	    lMt = length(indmaxtmp);
	    if lmt > 0
	      indminps(k,1:lmt,i) = indmintmp + max([(limpsl(k,i)-1),1])-1;
	    end
	    if lMt > 0
	      indmaxps(k,1:lMt,i) = indmaxtmp + max([(limpsl(k,i)-1),1])-1;
	    end
	    if lmt < size(indminps,2)
	      indminps(k,length(indmintmp)+1:end,i) = 0;
	    end
	    if lMt < size(indmaxps,2)
	      indmaxps(k,length(indmaxtmp)+1:end,i) = 0;
	    end
	    needextrps(k,i) = 0;
	    
	  end
	  
	  if stoprps(k,i) >= LX
	    lafinps(k,i) = 1;
	    needextrps(k,i) = 0;
	  end
	end
	
	if ~suspps(k,i)
	  curindminps = indminps(k,find(indminps(k,:,i) >= limpsl(k,i)),i);
	  curindmaxps = indmaxps(k,find(indmaxps(k,:,i) >= limpsl(k,i)),i);
	  nemps = length(curindminps) + length(curindmaxps);
	  
	end
      
	% loop of block (pre)sifting
	while (~needextrps(k,i) | lafinps(k,i)) & ~stoptestps(k,i) & ~waittest & ~suspps(k,i)

	  if nemps < 3 & lafinps(k,i)
	    stoptestps(k,:) = 1;
	    stoptest(k) = 1;
	    m(k,:) = mps(k,:,i);
	    if i > 1
	      m(k+1,:) = x - sum(m(1:k,:));
	    end
	    break
	  end
	  
	  if limpsl(k,i) == 1
	    startps(k,i) = 1;
	  else
	    startps(k,i) = stopps(k,i);
	  end
	  if lafinps(k,i)
	    stopps(k,i) = LX;
	    stoptestps(k,i) = 1;
	  else
	    stopps(k,i) = min(curindminps(max([1,end - LARGTRANSPS+1])),curindmaxps(max([1,end - LARGTRANSPS+1]))); % si ~lafinps(k,i)
	  end
	  
      if startps(k,i) == stopps(k,i)
          pause   
          needextrps(k,i) = 1;
          break
      end
      
	  lmr = sum(curindminps > stopps(k,i));
	  lMr = sum(curindmaxps > stopps(k,i));
	  nemrps(k,i) = lmr + lMr;
	  
	  if nemrps(k,i) < 8*LARGTRANSPS
	    needextrps(k,i) = 1;
	  end
	  
	  
	  if limpsl(k,i) == 1
	    margeml = 0;
	    margeMl = 0;
	      
       if curindmaxps(1) < curindminps(1)
		if mps(k,1,i) > mps(k,curindminps(1),i)
	      lmax = fliplr(curindmaxps(2:min(end,NBSYM+1)));
	      lmin = fliplr(curindminps(1:min(end,NBSYM)));
	      lsym = curindmaxps(1);
		else
	      lmax = fliplr(curindmaxps(1:min(end,NBSYM)));
	      lmin = [fliplr(curindminps(1:min(end,NBSYM-1))),1];
	      lsym = 1;
		end
	      else
		if mps(k,1,i) < mps(k,curindmaxps(1),i)
	      lmax = fliplr(curindmaxps(1:min(end,NBSYM)));
	      lmin = fliplr(curindminps(2:min(end,NBSYM+1)));
	      lsym = curindminps(1);
		else
	      lmax = [fliplr(curindmaxps(1:min(end,NBSYM-1))),1];
	      lmin = fliplr(curindminps(1:min(end,NBSYM)));
	      lsym = 1;
		end
	      end
	      
	    tlmin = 2*t(lsym)-t(lmin);
	    tlmax = 2*t(lsym)-t(lmax);
	    
    % in case symmetrized parts do not extend enough
    if tlmin(1) > t(1) | tlmax(1) > t(1)
      if lsym == curindmaxps(1)
	lmax = fliplr(curindmaxps(1:min(end,NBSYM)));
      else
	lmin = fliplr(curindminps(1:min(end,NBSYM)));
      end
      if lsym == 1
	error('bug')
      end
      lsym = 1;
      tlmin = 2*t(lsym)-t(lmin);
      tlmax = 2*t(lsym)-t(lmax);
    end
	    	    
	  else % limpsl ~= 1
	    
	    lmin = curindminps(find(curindminps <= startps(k,i)));
	    lmax = curindmaxps(find(curindmaxps <= startps(k,i)));
	    if length(lmin) < 5 |length(lmax) < 5
	      error('souci')
	    end
	    tlmin = t(lmin);
	    tlmax = t(lmax);
	    margeml = length(lmin);
	    margeMl = length(lmax);
	    
	  end % if limpsl...
	    
	  if lafinps(k,i)
	    
	    margemr = 0;
	    margeMr = 0;
	      
	      if curindmaxps(end) < curindminps(end)
		if mps(k,LX,i) > mps(k,curindmaxps(end),i)
	rmax = [LX,fliplr(curindmaxps(max(end-NBSYM+2,1):end))];
	rmin = fliplr(curindminps(max(end-NBSYM+1,1):end));
	rsym = LX;
		else
	rmax = fliplr(curindmaxps(max(end-NBSYM+1,1):end));
	rmin = fliplr(curindminps(max(end-NBSYM,1):end-1));
	rsym = curindminps(end);
		end
	      else
		if mps(k,LX,i) < mps(k,curindminps(end),i)
 	rmax = fliplr(curindmaxps(max(end-NBSYM+1,1):end));
	rmin = [LX,fliplr(curindminps(max(end-NBSYM+2,1):end))];
	rsym = LX;
		else         
	rmax = fliplr(curindmaxps(max(end-NBSYM,1):end-1));
	rmin = fliplr(curindminps(max(end-NBSYM+1,1):end));
	rsym = curindmaxps(end);
		end
	      end
	      
	    trmin = 2*t(rsym)-t(rmin);
	    trmax = 2*t(rsym)-t(rmax);

    if trmin(end) < t(LX) | trmax(end) < t(LX)
      if rsym == indmax(end)
	rmax = fliplr(curindmaxps(max(end-NBSYM+1,1):end));
      else
	rmin = fliplr(curindminps(max(end-NBSYM+1,1):end));
      end
      if rsym == LX
 	error('bug')
      end
      rsym = LX;
      trmin = 2*t(rsym)-t(rmin);
      trmax = 2*t(rsym)-t(rmax);
    end

        
	  else % lafin(k,i) ~= 1
	    
	    rmin = curindminps(end - LARGTRANSPS + 1:end);
	    rmax = curindmaxps(end - LARGTRANSPS + 1:end);
	    trmin = t(rmin);
	    trmax = t(rmax);
	    margemr = length(rmin);
	    margeMr = length(rmax);
	    
    end
	  
	  mpslmax = mps(k,lmax,i);
	  mpslmin = mps(k,lmin,i);
	  mpsrmax = mps(k,rmax,i);
	  mpsrmin = mps(k,rmin,i);
	  
	  
	  envmax = interp1([tlmax t(curindmaxps((margeMl+1):(end-margeMr))) trmax],[mpslmax mps(k,curindmaxps((margeMl+1):(end-margeMr)),i) mpsrmax],t(startps(k,i):stopps(k,i)),'spline');
	  envmin = interp1([tlmin t(curindminps((margeml+1):(end-margemr))) trmin],[mpslmin mps(k,curindminps((margeml+1):(end-margemr)),i) mpsrmin],t(startps(k,i):stopps(k,i)),'spline');
	  
	  envmoy = (envmax + envmin)/2;
	  

	  
	  if i == NBPRESIFT
	    m(k,startps(k,i):stopps(k,i)) = mps(k,startps(k,i):stopps(k,i),i) - envmoy;
	    fin(k) = stopps(k,i);
	    susp(k) = 0;
	    if startps(k,i) == 1
	      nbmodes_psdone = k;
	    end
	  else
	    if nbstartedpresift(k) < i+1
	      nbstartedpresift(k) = i+1;
	    end
	    mps(k,startps(k,i):stopps(k,i),i+1) = mps(k,startps(k,i):stopps(k,i),i) - envmoy;
	    finps(k,i+1) = stopps(k,i);
	    suspps(k,i+1) = 0;
	  end
	  
	  limpsl(k,i) = min([curindminps(max(end - 15*LARGTRANSPS,1)),curindmaxps(max(end - 15*LARGTRANSPS,1))]);
	  
	  % display
	  if tst
	    figure(figures(3))
	    subplot(2*NBPRESIFT,1,2*i-1)
	    plot(t(stopps(k,i):finps(k,i)),mps(k,stopps(k,i):finps(k,i),i))
	    hold on
	    plot(t(1:startps(k,i)),mps(k,1:startps(k,i),i),'k')
	    plot(t(startps(k,i):stopps(k,i)),mps(k,startps(k,i):stopps(k,i),i),'r')
	    plot(t(startps(k,i):stopps(k,i)),envmax,'k--')
	    plot(t(startps(k,i):stopps(k,i)),envmin,'k--')
	    axis([1,LX,min(mps(k,1:finps(k,i),i)),max(mps(k,1:finps(k,i),i))])
	    hold off
	    title(['mode ',int2str(k),'. Presifting # ',int2str(i)])
	    subplot(2*NBPRESIFT,1,2*i)
	    if i == NBPRESIFT
	      plot(t(startps(k,i):finps(k,i)),m(k,startps(k,i):finps(k,i)))
	      hold on
	      plot(t(1:startps(k,i)),m(k,1:startps(k,i)),'k')
	      hold off
	      axis([1,LX,min(m(k,1:finps(k,i))),max(m(k,1:finps(k,i)))])
	    else
	      plot(t(startps(k,i):stopps(k,i)),mps(k,startps(k,i):stopps(k,i),i+1))
	      hold on
	      plot(t(1:startps(k,i)),mps(k,1:startps(k,i),i+1),'k')
	      plot(t(stopps(k,i):finps(k,i+1)),mps(k,stopps(k,i):finps(k,i+1),i+1))
	      hold off
	      axis([1,LX,min(mps(k,1:finps(k,i+1),i+1)),max(mps(k,1:finps(k,i+1),i+1))])
	    end
	    pause(0.01)
	  end
	    
	    
	  
	  
	  
	end % petit while
      end % for
    end % grand while 
    
    
    % interrupts loop on a mode for processing to the next one    
    waittest = 0;
    
    
    
    while ~stoptest(k) & ~waittest & ~susp(k)	

	
      if needextr(k)
	if fin(k) > 0
	  stopr(k) = fin(k);
	  [indm,indM] = extr(m(k,max([(startl(k)-1),1]):stopr(k)));
	else
	  indm = [];
	  indM  = [];
	  susp(k) = 1;
	end
    
	
	if sum([(indm > stop(k)),(indM > stop(k))]) < 4*LARGTRANS & fin(k) < LX %stopr(k) >= fin(k) & fin(k) < LX
	  susp(k) = 1;
	else
	  startr(k) = stopr(k) + 1;
	  [indmintmp,indmaxtmp] = extr(m(k,max([(startl(k)-1),1]):stopr(k)));
	  lmt = length(indmintmp);
	  lMt = length(indmaxtmp);
	  indmin(k,1:lmt) = indmintmp + max([(startl(k)-1),1])-1;
	  indmax(k,1:lMt) = indmaxtmp + max([(startl(k)-1),1])-1;
	  if lmt < size(indmin,2)
	    indmin(k,length(indmintmp)+1:end) = 0;
	  end
	  if lMt < size(indmax,2)
	    indmax(k,length(indmaxtmp)+1:end) = 0;
	  end
	  needextr(k) = 0;
	  if k==2
	  end
	end
	
	if stopr(k) >= LX
	  lafin(k) = 1;
	  stop(k) = LX;
	  needextr(k) = 0;
	end
      end
      
      if ~susp(k)
	curindmin = indmin(k,find(indmin(k,:) >= startl(k)));
	curindmax = indmax(k,find(indmax(k,:) >= startl(k)));
	nem = length(curindmin) + length(curindmax);
	
	lml(k) = sum(curindmin < start(k));
	lMl(k) = sum(curindmax < start(k));
      end

      % loop of local on-line sifting
      while (~needextr(k) | lafin(k)) & ~stoptest(k) & ~waittest & ~susp(k)
		
	if nem < 3 & lafin(k) == 1
	  stoptest(k) = 1;
	  if nbit(k) > 1 
	    m(k+1,:) = x - sum(m(1:k,:));
	  end
	  break
	end
	
	dm = diff(curindmin);
	dM = diff(curindmax);
	
	
	if lml(k)+lMl(k) < 4 | (start(k) < max([curindmin(min(LARGTRANS,end)),curindmax(min(LARGTRANS,end))]) + mean([dm(1:min([end,5])),dM(1:min([end,5]))]))
	  margeml = 0;
	  margeMl = 0;
	    
	    if curindmax(1) < curindmin(1)
	      if m(k,1) > m(k,curindmin(1))
	lmax = fliplr(curindmax(2:min(end,NBSYM+1)));
	lmin = fliplr(curindmin(1:min(end,NBSYM)));
	lsym = curindmax(1);
	      else
	lmax = fliplr(curindmax(1:min(end,NBSYM)));
	lmin = [fliplr(curindmin(1:min(end,NBSYM-1))),1];
	lsym = 1;
	      end
	    else
	      if m(k,1) < m(k,curindmax(1))
	lmax = fliplr(curindmax(1:min(end,NBSYM)));
	lmin = fliplr(curindmin(2:min(end,NBSYM+1)));
	lsym = curindmin(1);
	      else
	lmax = [fliplr(curindmax(1:min(end,NBSYM-1))),1];
	lmin = fliplr(curindmin(1:min(end,NBSYM)));
	lsym = 1;
	      end
	    end
	    
	  tlmin = 2*t(lsym)-t(lmin);
	  tlmax = 2*t(lsym)-t(lmax);

      % in case symmetrized parts do not extend enough
    if tlmin(1) > t(1) | tlmax(1) > t(1)
      if lsym == curindmax(1)
	lmax = fliplr(curindmax(1:min(end,NBSYM)));
      else
	lmin = fliplr(curindmin(1:min(end,NBSYM)));
      end
      if lsym == 1
	error('bug')
      end
      lsym = 1;
      tlmin = 2*t(lsym)-t(lmin);
      tlmax = 2*t(lsym)-t(lmax);
    end

      debfen = startl(k);
	  
	else
	  lmin = curindmin(1:LARGTRANS);
	  margeml = LARGTRANS;
	  lmax = curindmax(1:LARGTRANS);
	  margeMl = LARGTRANS;
	  debfen = max([lmin,lmax]);
	  tlmin = t(lmin);
	  tlmax = t(lmax);
	end
	
	if lafin(k) == 1
	  margemr = 0;
	  margeMr = 0;

        if curindmax(end) < curindmin(end)
	      if m(k,LX) > m(k,curindmax(end))
	rmax = [LX,fliplr(curindmax(max(end-NBSYM+2,1):end))];
	rmin = fliplr(curindmin(max(end-NBSYM+1,1):end));
	rsym = LX;
	      else
	rmax = fliplr(curindmax(max(end-NBSYM+1,1):end));
	rmin = fliplr(curindmin(max(end-NBSYM,1):end-1));
	rsym = curindmin(end);
	      end
	    else
	      if m(k,LX) < m(k,curindmin(end))
	rmax = fliplr(curindmax(max(end-NBSYM+1,1):end));
	rmin = [LX,fliplr(curindmin(max(end-NBSYM+2,1):end))];
	rsym = LX;
	      else         
	rmax = fliplr(curindmax(max(end-NBSYM,1):end-1));
	rmin = fliplr(curindmin(max(end-NBSYM+1,1):end));
	rsym = curindmax(end);
	      end
	    end
	  trmin = 2*t(rsym)-t(rmin);
	  trmax = 2*t(rsym)-t(rmax);

    if trmin(end) < t(LX) | trmax(end) < t(LX)
      if rsym == curindmax(end)
	rmax = fliplr(curindmax(max(end-NBSYM+1,1):end));
      else
	rmin = fliplr(curindmin(max(end-NBSYM+1,1):end));
      end
      if rsym == LX
 	error('bug')
      end
      rsym = LX;
      trmin = 2*t(rsym)-t(rmin);
      trmax = 2*t(rsym)-t(rmax);
    end

    
      finfen(k) = LX;
	  
	else
	  rmin = curindmin((end - LARGTRANS + 1):end);
	  rmax = curindmax((end - LARGTRANS + 1):end);
	  margemr = LARGTRANS;
	  margeMr = LARGTRANS;
	  finfen(k) = min([rmin,rmax]);
	  trmin = t(rmin);
	  trmax = t(rmax);
	end
	    
	mlmin = m(k,lmin);
	mlmax = m(k,lmax);
	mrmin = m(k,rmin);
	mrmax = m(k,rmax);
    
	lcur = finfen(k) - debfen + 1;
	
	
	envmax = interp1([tlmax t(curindmax(margeMl+1:(end-margeMr))) trmax],[mlmax m(k,curindmax(margeMl+1:(end-margeMr))) mrmax],t(debfen:finfen(k)),'spline');
	envmin = interp1([tlmin t(curindmin(margeml+1:(end-margemr))) trmin],[mlmin m(k,curindmin(margeml+1:(end-margemr))) mrmin],t(debfen:finfen(k)),'spline');
		      
	envmoy = (envmax + envmin)/2;

	
	amp = abs(envmax-envmin)/2;
	mamp = max(amp);
	sx = abs(envmoy)./amp;
	s = mean(sx);
	
	
	
    % definition of f equal to 1 where sifting is needed and decays 
    % fast to 0 in the neighbourhood

	badsx = (sx > sd & amp > mamp/tol) | (sx > sd2 & amp > mamp/tol2);
	d = diff([0 badsx 0]);
	debs = find(d==1);
	fins = find(d==-1)-1;
	if length(debs)~=length(fins)
	  error('pb avec les composantes connexes')
	end
	indextr = sort([curindmin curindmax]);
	d = diff(indextr);
	connexe = [];
	lc = length(debs);
	f = [];
	if lc > 0
	  connexe(lc,lcur) = 0;
	  for i = 1:lc
	    connexe(i,debs(i):fins(i)) = 1;
	    indp = min([lc-1,length(find(indextr < debs(i)))]);
	    inds = max([2,length(find(indextr <= fins(i)))+1]);
	    llarg = mean(d(min([max([nem-1,1]),max([1,indp])]):max([1,min([nem-1,indp + 2])])));
	    rlarg = mean(d(min([max([nem-1,1]),max([1,inds - 1])]):max([1,min([nem-1,inds + 1])])));
	    larg = mean([rlarg,llarg]);
	    if indp == 0
	      if inds ~= (nem+1)
		larg = (indextr(inds)-indextr(inds-1));
	      else
		larg = round(lcur/2);
	      end
	    else
	      if  inds ~= (nem+1)
		larg = max([(indextr(inds)-indextr(inds-1)),(indextr(indp+1)-indextr(indp))]);
	      else
		larg =(indextr(indp+1)-indextr(indp));
	      end
	    end
	    larg = 2 * max([round((fins(i)-debs(i))/4),larg,LARGMIN]);
	    w(1:round(larg/2)) = [1:round(larg/2)]/round(larg/2);
	    w((2*larg+2-round(larg/2)):(2*larg+1)) = fliplr(w(1:round(larg/2)));
	    w(round(larg/2):(2*larg+2-round(larg/2))) = 1;
	    indd = max(1,debs(i)-larg);
	    indf = min(lcur,fins(i)+larg);
	    connexe(i,indd:debs(i)) = w((larg+1-debs(i)+indd):(larg+1));
	    connexe(i,fins(i):indf) = w((larg+1):(larg+1+indf-fins(i)));
	  end
	  
	  f = max(connexe,[],1);
	  
	else
	  f(lcur) = 0;
	end
	    
    %                                   ___
    % definition of f2, window of form /   \ used as a multiplier on f
	f2 = [];
	f2((start(k)-debfen+1):(stop(k)-debfen+1)) = 1;
	if debfen == 1
	  dl = round(max([mean([dm(1:min([end,3])),dM(1:min([end,3]))]),start(k)/2]));
	  f2(1:start(k)) = max([dl-start(k)+1:dl]/dl,0);
	else
	  f2(1:(start(k)-debfen+1)) = [0:(start(k)-debfen)]/(start(k)-debfen);
	end
	if lafin(k) == 1
	  f2(stop(k)-debfen+1:LX-debfen+1) = 1;
	else
	  f2((stop(k)-debfen+1):(finfen(k)-debfen+1)) = fliplr([0:(finfen(k)-stop(k))])/(finfen(k)-stop(k));
	end
	
	f1 = f;
	f = f.*f2;
	    
	    
    mp = m(k,:);
	m(k,debfen:finfen(k)) = m(k,debfen:finfen(k)) - f.*envmoy;
	    
	% display
	    
	if tst == 1
	  figure(figures(1))
	  subplot(4,1,1)
	  plot(t(1:fin(k)),mp(1:fin(k)));hold on;
	  axis([t(1),t(end),-max(abs(mp)),max(abs(mp))])
	  title(['IMF ',int2str(k),';   iteration ', int2str(nbit(k)),' before sifting']);
	  plot(t(debfen:finfen(k)),envmax,'--k');plot(t(debfen:finfen(k)),envmin,'--k');plot(t(debfen:finfen(k)),envmoy,'r');
	  hold  off
      set(gca,'XTick',[])
	  subplot(4,1,2)
	  plot(t(debfen:finfen(k)),sx);
      hold on
      plot(t,sdt,'--r')
      plot(t,sd2t,':k')
      title('stop parameter')
      set(gca,'XTick',[])
	  axis([t(1),t(end),0,max(sx)])
	  hold off
	  subplot(4,1,3);
	  plot(t(debfen:finfen(k)),f);
	  hold on;
	  plot(t(debfen:finfen(k)),f2,'--k');
	  hold off
	  axis([t(1),t(end),0,1.1])
	  hold off
      title('window')
      subplot(4,1,4);
	  plot(t,m(k,:));
	  axis([t(1),t(end),-max(abs(m(k,:))),max(abs(m(k,:)))])
	  title(['IMF ',int2str(k),';   iteration ', int2str(nbit(k)),' after sifting']);
      set(gca,'XTick',[])
      pause(0.01)
	  % clf
	end
	
	
	
	% if start change
	if sum(badsx(round((start(k)-debfen+1)/2):min([(start(k)+5*PAS-debfen+1),(stop(k) - debfen +1)])))==0
	  if start(k) >= LX
	    stoptest(k) = 1;
	    startl(k) = LX;
	    suspps(k+1,1) = 0;
	    disp(['mode ',int2str(k),' terminÈ'])
	    if nbmodes >= k + 1
	      finpsp = finps(k+1,1);
	    else
	      finpsp = 0;
	      nbmodes = nbmodes +1;
	      startl(k+1) = 1;
	      start(k+1) = 1;
	      stop(k+1) = min(PAS+1,LX);
	      startr(k+1) = 1;
	      stopr(k+1) = 1;
	      lafin(k+1) = 0;
	      stoptest(k+1) = 0;
	      susp(k+1) = 0;
	      needextrps(k+1,1:NBPRESIFT) = 1;
	      needextr(k+1) = 1;
	      nbit(k+1) = 0;
	      lm(k+1) = 0;
	      lM(k+1) = 0;
	      nemr(k+1) = 0;
	      finfen(k+1) = 1;
	      fin(k+1) = 0;
	      limpsl(k+1,:) = 1;
	      startps(k+1,:) = 1;
	      stopps(k+1,:) = 1;
	      lafinps(k+1,NBPRESIFT) = 0;
	      suspps(k+1,NBPRESIFT) = 0;
	      stoptestps(k+1,NBPRESIFT) = 0;
	      nbstartedpresift(k+1) = 1;
	      m(k+1,:) = 0;
	    end
	    finps(k+1,1) = LX;
	    mps(k+1,finpsp+1:finps(k+1,1),1) = x(finpsp+1:finps(k+1,1))-sum(m(1:k,finpsp+1:finps(k+1,1)),1);
	  end
	  tmp = find(f(start(k)-debfen+1:stop(k)-debfen+1) > 0);
	  if length(tmp) == 0
	    if lafin(k)
	      start(k) = LX;
	    else
	      start(k) = min(round((stop(k) + start(k))/2),LX);
	    end
	  else
	    start(k) = min(start(k)+round((tmp(1)-1)/2),LX);
	  end
	  
	  lmlp = lml(k);
	  lMlp = lMl(k);
	  
	  lml(k) = sum(curindmin < start(k));
	  lMl(k) = sum(curindmax < start(k));
	  if lml(k)+lMl(k) > 4*LARGTRANS + 1
	    tmp = lml(k)+lMl(k) - 4*LARGTRANS;
	    startlp = startl(k);
	    if ~stoptest(k)
	      if curindmin(1) > curindmax(1)
		startl(k) = curindmax(ceil(tmp/2));
	      else
		startl(k) = curindmin(ceil(tmp/2));
	      end
	    else
	      startl(k) = LX;
	    end
	    if startlp == 1 & nbmodes < k + 1
	      nbmodes = nbmodes +1;
	      startl(k+1) = 1;
	      start(k+1) = 1;
	      stop(k+1) = min(PAS+1,LX);
	      startr(k+1) = 1;
	      stopr(k+1) = 1;
	      lafin(k+1) = 0;
	      stoptest(k+1) = 0;
	      susp(k+1) = 0;
	      needextrps(k+1,1:NBPRESIFT) = 1;
	      needextr(k+1) = 1;
	      nbit(k+1) = 0;
	      lm(k+1) = 0;
	      lM(k+1) = 0;
	      nemr(k+1) = 0;
	      fin(k+1) = 0;
	      finfen(k+1) = 1;
	      limpsl(k+1,:) = 1;
	      startps(k+1,:) = 1;
	      stopps(k+1,:) = 1;
	      finps(k+1,:) = 0;
	      lafinps(k+1,NBPRESIFT) = 0;
	      suspps(k+1,NBPRESIFT) = 0;
	      stoptestps(k+1,NBPRESIFT) = 0;
	      nbstartedpresift(k+1) = 1; 
	      m(k+1,:) = 0;
	    end
	    
	    finpsp = finps(k+1,1);
	    finps(k+1,1) = startl(k);
	    suspps(k+1) = 0;
	    mps(k+1,finpsp+1:finps(k+1,1),1) = x(finpsp+1:finps(k+1,1)) - sum(m(1:k,finpsp+1:finps(k+1,1)),1);
	    curindmin = indmin(k,find(indmin(k,:) >= startl(k)));
	    curindmax = indmax(k,find(indmax(k,:) >= startl(k)));
	  end
	end      
	
	nsteps = nsteps+1;
	nbit(k) = nbit(k) + 1;
	if nsteps == STEP
	  waittest = 1;
	end
	    
	    

	stop(k) = min([round((4*stop(k) + finfen(k))/5),fin(k)]);

	
	[indmintmp,indmaxtmp] = extr(m(k,max([(startl(k)-1),1]):stopr(k)));
	lmt = length(indmintmp);
	lMt = length(indmaxtmp);
	indmin(k,1:lmt) = indmintmp + max([(startl(k)-1),1])-1;
	indmax(k,1:lMt) = indmaxtmp + max([(startl(k)-1),1])-1;
	if lmt < size(indmin,2)
	  indmin(k,length(indmintmp)+1:end) = 0;
	end
	if lMt < size(indmax,2)
	  indmax(k,length(indmaxtmp)+1:end) = 0;
	end
	curindmin = indmin(k,find(indmin(k,:) >= startl(k)));
	curindmax = indmax(k,find(indmax(k,:) >= startl(k)));
	nem = length(curindmin) + length(curindmax);
	
	lmr = sum(curindmin > stop(k));
	lMr = sum(curindmax > stop(k));
	nemr(k) = lmr + lMr;
	
	
    
	if nemr(k) < (4*LARGTRANS) & ~lafin(k)
	  needextr(k) = 1;
	end
	
	    
      end % inner loop of on-line sifting

    end % loop of on-line sifting
	
    % display of mode processed by block sifting
    
    if tst2 == 1
      figure(figures(2))
      if k <= nbmodes_psdone
	subplot(nbmodes_psdone,1,k)
	plot(t(startl(k):finfen(k)),m(k,startl(k):finfen(k)),'r');hold on;
	plot(t(1:startl(k)),m(k,1:startl(k)),'k');%hold on;
	plot(t(finfen(k):fin(k)),m(k,finfen(k):fin(k)));hold off;
	axis([t(1),t(end),-max(abs(m(k,1:fin(k)))),max(abs(m(k,1:fin(k))))])
      
      end
      pause(0.01)
    end
    
    
  end % for i = 1:nbmodes ...
end % global loop
    
    
imf = m;

ort = io(x,imf);

if tst | tst2
  disp('press any key to continue')
  pause
end

if tst
    close(figures(1))
    close(figures(3))
end
if tst2
    close(figures(2))
end
