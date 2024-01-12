function [ output_args ] = setGcf( in )
%SET_GCF Summary of this function goes here
%   Detailed explanation goes here


if isstr(in) || (length(in) == 3)
    set(gcf, 'color', in)
else
    set(gcf, 'color', in * [1,1,1])
end

