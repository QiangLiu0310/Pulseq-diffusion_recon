function [ output_args ] = s( input_args, dim_choose )
%S Summary of this function goes here
%   Detailed explanation goes here

if nargin > 1
    output_args = size(input_args, dim_choose);
else
    output_args = size(input_args);
end

end

