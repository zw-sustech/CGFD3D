% INDEX_EMD.M list of functions in the EMD package
% 
% type help function_name for more information on a specific function
% 
% Empirical Mode Decomposition
%
%   emd          - computes EMD and bivariate/complex EMD with various options
%   emd_local    - computes local EMD variation
%   emd_online   - computes on-line EMD variation. Note that it does not truly
%                  apply on-line: the function is only a demonstration.
%   emdc         - fast implementation for EMD with Cauchy-like stopping criterion
%                  (requires compilation, see make_emdc function)
%   emdc_fix     - fast implementation for EMD with predefined number of iterations
%                  (requires compilation, see make_emdc function)
%   cemdc        - fast implementation for bivariate/complex EMD (first algorithm)
%                  with Cauchy-like stopping criterion (requires compilation,
%                  see make_emdc function)
%   cemdc_fix    - fast implementation for bivariate/complex EMD (first algorithm)
%                  with predefined number of iterations (requires compilation,
%                  see make_emdc function)
%   cemdc2       - fast implementation for bivariate/complex EMD (second algorithm)
%                  with Cauchy-like stopping criterion (requires compilation,
%                  see make_emdc function)
%   cemdc2_fix   - fast implementation for bivariate/complex EMD (second algorithm)
%                  with predefined number of iterations (requires compilation,
%                  see make_emdc function)
% 
% Utilities
% 
%   install_emd   - setup Matlab's path and compile the C codes.
%   uninstall_emd - revert the modifications made by install_emd and remove the 
%                   files (optional).
%   make_emdc     - compile all C codes
%   emd_visu      - visualization of EMD
%   cemd_visu     - visualization of bivariate/complex EMD (automatically called
%                   by emd_visu when the input is complex)
%   cenvelope     - compute envelope curves for bivariate/complex EMD
%   cemd_disp     - visualization of envelope curves and tube envelope
%   plot3c        - plot a complex vector in 3 dimensions
%   plotc         - plot the projection of a complex vector on a variable direction
%   dirstretch    - directional stretching of a complex vector
%   hhspectrum    - compute Hilbert-Huang spectrum (need the Time-Frequency Toolbox
%                   http://tftb.nongnu.org)
%   toimage       - transform a spectrum made of 1D functions (e.g., output of
%                   "hhspectrum") in an 2D image
%   disp_hhs      - display the image output of "toimage" as a Hilbert-Huang spectrum
%   addtag        - add a tag to a graphic object (uses the Tag property as a list
%                   of keywords or "tags")
%   rmtag         - remove a tag from a graphic object (uses the Tag property as
%                   a list of keywords or "tags")
%   hastag        - test whether a graphic object has a specific tag (uses the Tag
%                   property as a list of keywords or "tags")
%   findtag       - find objects having a specific tag (uses the Tag property as
%                   a list of keywords or "tags")
%   
% Examples from G. Rilling, P. Flandrin and P. Gonçalves,
%   "On Empirical Mode Decomposition and its algorithms"
%   IEEE-EURASIP Workshop on Nonlinear Signal and Image Processing
%   NSIP-03, Grado (I), June 2003
%   
%   emd_fmsin         - Fig. 1: a 3-component example (need the Time-Frequency 
%                       Toolbox http://tftb.nongnu.org)
%   emd_triang        - Fig. 2: another 3-component example
%   emd_sampling      - Fig. 3: effect of sampling on 1 tone
%   emd_separation    - Fig. 4: separation of 2 tones
%   ex_online         - Sect 3.4: the way emd_online.m works
%   triangular_signal - subroutine called by emd_triang (formerly triang.m)
%   
% Examples from G. Rilling, P. Flandrin, P. Gonçalves and J. M. Lilly,
%   "Bivariate Empirical Mode Decomposition",
%   Signal Processing Letters (submitted)
% 
%   bivariate_EMD_principle        - Fig. 1: principle of the bivariate/complex EMD
%   bivariate_EMD_mean_definitions - Fig. 2: definition of the mean for each algorithm. 
%                                    Also allows to test other signals and parameter sets.
%   bivariate_EMD_illustration     - Fig. 3: illustration of the bivariate EMD
%                                    on an oceanographic float position record


help index_emd
  