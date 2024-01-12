function [ res ] = apodize3d( img, apodization_para )
%APODIZE Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
    apodization_para=0.2;
end

kcomb=img;
kapodize = mrir_filter_raw_apodize_1d( mrir_filter_raw_apodize_1d(kcomb, 1, apodization_para),  2, apodization_para);
res=kapodize;
end

