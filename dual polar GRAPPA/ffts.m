function y = ffts(x,dim)
% FFTS Shift zero-frequency component to center of spectrum.
%   For vectors, FFTSHIFT(X) swaps the left and right halves of
%   X.  For matrices, FFTSHIFT(X) swaps the first and third
%   quadrants and the second and fourth quadrants.  For N-D
%   arrays, FFTSHIFT(X) swaps "half-spaces" of X along each
%   dimension.
%
%   FFTSHIFT(X,DIM) applies the FFTSHIFT operation along the 
%   dimensions specified in the array DIM.
%
%   FFTSHIFT is useful for visualizing the Fourier transform with
%   the zero-frequency component in the middle of the spectrum.
%
%   Class support for input X:
%      float: double, single
%
%   See also IFFTSHIFT, FFT, FFT2, FFTN, CIRCSHIFT.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 5.11.4.3 $  $Date: 2004/03/09 16:16:18 $

numDims = ndims(x);  
if nargin > 1
  
else,
  dim = 1:numDims;
end;

idx = cell(1, numDims);
for k = 1:numDims
  if ( sum(k == dim) > 0 ),
    m = size(x, k);
    p = ceil(m/2);
    idx{k} = [p+1:m 1:p];
  else
    idx{k} = 1:size(x,k);
  end;
end


% Use comma-separated list syntax for N-D indexing.
y = x(idx{:});
