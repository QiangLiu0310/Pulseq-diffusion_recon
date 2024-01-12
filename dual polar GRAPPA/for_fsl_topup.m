clc;clear all;close all;
% load('pos_b0.mat'); load('pos_b1.mat'); load('pos_b2.mat'); load('pos_b3.mat')
% load('neg_b0.mat');load('neg_b1.mat');load('neg_b2.mat');load('neg_b3.mat')
% dir=4;
% pos=zeros(size(posb0,1),size(posb0,2),size(posb0,3),dir);
% neg=zeros(size(posb0,1),size(posb0,2),size(posb0,3),dir);
% for i=1:dir
%     pos(:,:,:,i)=eval(strcat('posb',num2str(i-1)));
%     neg(:,:,:,i)=eval(strcat('negb',num2str(i-1)));
% end   
% load('blipup.mat')
% load('blipdown.mat')
% load('D:\dpg_buda_postpro\dpg_buda_1027\blipum.mat')
 load('blip_down.mat')
 load('blip_up.mat')
blipd=squeeze(sos(fft2c(blip_down(:,:,:,2)),3));
blipu=squeeze(sos(fft2c(blip_up(:,:,:,2)),3));
blipd=flipdim(permute(repmat(blipd,1,1,1,4),[2 1 4 3]),1);
blipu=flipdim(flipdim(permute(repmat(blipu,1,1,1,4),[2 1 4 3]),1),2);

% I=blipu(:,:,1);
% figure,imshow(I,[])
% h=imfreehand;
% mask1=createMask(h);
% 
%   I=blipu(:,:,1);
% figure,imshow(I,[])
% h=imfreehand;
% mask2=createMask(h);
%  mask=double(mask1)+double(mask2);
%   blipu=blipu.*mask;
%  
%  save('blipum.mat','mask')


blipdnii= make_nii(blipd(:,:,:,:));
save_nii(blipdnii,'blip_down.nii')
% gzip('blip_down.nii')
blipunii= make_nii(blipu(:,:,:,:));
save_nii(blipunii,'blip_up.nii')
% gzip('blip_up.nii')

files = gunzip('my_hifi_images.nii.gz')
nii=load_nii('my_hifi_images.nii')
A=nii.img;
files=gunzip('my_topup_results_fieldcoef.nii.gz')
nii=load_nii('my_topup_results_fieldcoef.nii')
% imshow3(A(:,:,1))
% A=imresize(A,[128 96]);
A=flipdim(A,1);
files=gunzip('both_b0.nii.gz')
nii=load_nii('both_b0.nii')
A=nii.img;