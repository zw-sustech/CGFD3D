%IO  computes the index of orthogonality
%
% ORT = IO(X,IMF)
%
%                       ________
%         _  |IMF(i,:).*IMF(j,:)|
%   ORT = \ _____________________
%         /
%         ¯        || X ||²
%        i~=j
%
% inputs : - X    : analyzed signal
%          - IMF  : empirical mode decomposition
function ort = io(x,imf)

lx = size(imf,2);
n = size(imf,1);

s = 0;

for i = 1:n
  for j =1:n
    if i~=j
      s = s + abs(sum(imf(i,:).*conj(imf(j,:)))/sum(x.^2));
    end
  end
end

ort = 0.5*s;
