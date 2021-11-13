%EMDC  computes Empirical Mode Decomposition
%
%
%   Syntax
%
%
% [IMF,NB_ITERATIONS]=EMDC(T,X,STOP_PARAMETERS,MAX_IMFS);
%
%
%   Description
%
%
% computes EMD according to [1] with stopping criterion for sifting in [2]:
%
%   mean of boolean array {(mean_amplitude)/(envelope_amplitude) > THRESHOLD} < TOLERANCE
%   &
%   |#zeros-#extrema|<=1
%
% inputs:	
%       - T: sampling times. If T=[], the signal is assumed uniformly sampled.
%       - X: analyzed signal
%       - STOP_PARAMETERS: parameters for the stopping criterion: 
%         if scalar the value is used to specify THRESHOLD only.
%         otherwise the vector should be: [THRESHOLD,TOLERANCE].
%         if STOP_PARAMETERS is unspecified or empty, default values are used: [0.05,0.05]
%       - MAX_IMFS: maximum number of IMFs to be extracted. If MAX_IMFS is
%         zero, empty or unspecified, the default behavior is to extract as
%         many IMFs as possible.
%         
% outputs: 
%		- IMF: intrinsic mode functions (IMFs) (last line = residual)
%		- NB_ITERATIONS: effective number of sifting iterations for each mode
%
%
%   Examples
%
%
% workspace: 
%  T: 1xN time instants
%  X: 1xN signal data 
%
%>>IMF = EMDC(T,X);
%>>[IMF,NB_IT] = EMDC([],X);
%>>IMF = EMDC(T,X,0.1);
%>>IMF = EMDC(T,X,[0.1,0.1]);
%>>[IMF,NB_IT] = EMDC([],X,[],4);
%
%
%   References
%
%
% [1] N. E. Huang et al., "The empirical mode decomposition and the
% Hilbert spectrum for non-linear and non stationary time series analysis",
% Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998
%
% [2] G. Rilling, P. Flandrin and P. Gonçalves
% "On Empirical Mode Decomposition and its algorithms",
% IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
% NSIP-03, Grado (I), June 2003
%
%
% See also
%  emd_visu (visualization),
%  emd (slow but has many options),
%  emdc_fix (fast implementation of EMD with different stopping criterion)
%  hhspectrum (Hilbert-Huang spectrum)
%
%
% G. Rilling, last modification: 3.2007
% gabriel.rilling@ens-lyon.fr
%
% code based on a student project by T. Boustane and G. Quellec, 11.03.2004
% supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
% email : pchainai@isima.fr).
