function matrixOut = smooth2(matrixIn, Nr, Nc)

% SMOOTH2.M: Smooths matrix data.

%			MATRIXOUT=SMOOTH2(MATRIXIN,Nr,Nc) smooths the data in MATRIXIN 
%           using a running mean over 2*N+1 successive points, N points on 
%           each side of the current point.  At the ends of the series 
%           skewed or one-sided means are used.  
%
%           Inputs: matrixIn - original matrix
%                   Nr - number of points used to smooth rows
%                   Nc - number of points to smooth columns
%           Outputs:matrixOut - smoothed version of original matrix
%
%           Remark: By default, if Nc is omitted, Nc = Nr.
%
%           Written by Kelly Hilands, October 2004
%           Applied Research Laboratory
%           Penn State University
%
%           Developed from code written by Olof Liungman, 1997
%			Dept. of Oceanography, Earth Sciences Centre
%			Göteborg University, Sweden
%			E-mail: olof.liungman@oce.gu.se


%Initial error statements and definitions
if nargin<2, error('Not enough input arguments!'), end

N(1) = Nr; 
if nargin<3 
    N(2) = N(1); 
else
    N(2) = Nc;
end

if length(N(1))~=1, error('Nr must be a scalar!'), end
if length(N(2))~=1, error('Nc must be a scalar!'), end

for dim = 1:2       %dim = 1, smooth rows; dim = 2, smooth cols
    
    %Transpose matrix to smooth columns (instead of rows)
    if dim == 2
        matrixIn = matrixOut';
    end

    [row,col] = size(matrixIn);

    %Initialize temporary vectors
    vectorOut = zeros(1,col);
    matrixOut = zeros(row,col);
    temp = zeros(2*N(dim)+1,col-2*N(dim));
    temp(N(dim)+1,:) = matrixIn(N(dim)+1:col-N(dim));

    %Smooth each vector in the matrix
    for i = 1:row
        vectorIn = matrixIn(i,:);
        
        for j = 1:N(dim)
          vectorOut(j) = mean(vectorIn(1:j+N(dim)));
          vectorOut(col-j+1) = mean(vectorIn(col-j-N(dim):col));
          temp(j,:) = vectorIn(j:col-2*N(dim)+j-1);
          temp(N(dim)+j+1,:) = vectorIn(N(dim)+j+1:col-N(dim)+j);
        end
        
        vectorOut(N(dim)+1:col-N(dim)) = mean(temp);
        matrixOut(i,:) = vectorOut;
    end
end

%Transpose matrix back to original orientation
matrixOut = matrixOut';




