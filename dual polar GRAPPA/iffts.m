function y = iffts(x,dim)
%IFFTSHIFT Inverse FFT shift.
%   For vectors, IFFTSHIFT(X) swaps the left and right halves of
%   X.  For matrices, IFFTSHIFT(X) swaps the first and third
%   quadrants and the second and fourth quadrants.  For N-D
%   arrays, IFFTSHIFT(X) swaps "half-spaces" of X along each
%   dimension.
%
%   IFFTSHIFT(X,DIM) applies the IFFTSHIFT operation along the 
%   dimension DIM.
%
%   IFFTSHIFT undoes the effects of FFTSHIFT.
%
%   Class support for input X: 
%      float: double, single
%
%   See also FFTSHIFT, FFT, FFT2, FFTN.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.6.4.3 $  $Date: 2004/03/09 16:16:25 $


numDims = ndims(x);
if nargin > 1

else
  dim = 1:numDims;
end;

idx = cell(1, numDims);
for k = 1:numDims
  if ( sum(k == dim) > 0 ),
    m = size(x, k);
    p = floor(m/2);
    idx{k} = [p+1:m 1:p];
  else,
    idx{k} = 1:size(x,k);
  end;
end


% Use comma-separated list syntax for N-D indexing.
y = x(idx{:});
