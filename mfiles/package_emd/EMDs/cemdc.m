%CEMDC2  bivariate Empirical Mode Decomposition, first algorithm
%
%
%   Syntax
%
%
% [IMF,NB_ITERATIONS]=CEMDC(T,X,STOP_PARAMETERS,MAX_IMFS,NDIRS);
%
%
%   Description
%
%
% computes bivariate EMD, first algorithm [1] with stopping criterion for
% sifting similar to the one proposed in [2]:
%
%   mean of boolean array {(mean_amplitude)/(envelope_amplitude) > THRESHOLD} < TOLERANCE
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
%       - NDIRS: number of directions used to compute the local mean.
%         If unspecified, the default value is 4.
%         rem: the actual number of directions (according to [1]) is 2*NDIRS
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
%>>IMF = CEMDC(T,X);
%>>[IMF,NB_IT] = CEMDC([],X);
%>>IMF = CEMDC(T,X,0.1);
%>>IMF = CEMDC(T,X,[0.1,0.1]);
%>>[IMF,NB_IT] = CEMDC([],X,[],4);
%
%
%   References
%
%
% [1] G. Rilling, P. Flandrin, P. Gonçalves and J. M. Lilly.,
% "Bivariate Empirical Mode Decomposition",
% Signal Processing Letters (submitted)
%
% [2] G. Rilling, P. Flandrin and P. Gonçalves
% "On Empirical Mode Decomposition and its algorithms",
% IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
% NSIP-03, Grado (I), June 2003
%
%
% See also
%(c)emd_visu (visualization),
%emd (slow but has many options),
%cemdc2_fix, cemdc, cemdc_fix (other fast implementations of bivariate EMD)
%
%
% G. Rilling, last modification: 3.2007
% gabriel.rilling@ens-lyon.fr
%
% code based on a student project by T. Boustane and G. Quellec, 11.03.2004
% supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
% email : pchainai@isima.fr).
