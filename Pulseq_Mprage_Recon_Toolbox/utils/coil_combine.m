function [ res ] = coil_combine( img, sens, chan_dim )
%COIL_COMBINE Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    chan_dim = ndims(img);
end


res = sum(conj(sens) .* img, chan_dim) ./ (eps + sum(abs(sens).^2, chan_dim));


end

