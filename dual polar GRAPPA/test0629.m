clc;clear all; close all;
% read in the raw EPI data (R=2)
load('shot1.mat')
k=shot1(:,:,:,1,1);
load('noblip.mat')

% compute the standard Nyquist ghost correction coefficients:
s=prof(1:3,:,:,1,1);
%s=permute(s,[2 1 3]);

% s= ifi(tmult(fif(s),diag(exp(j*4*pi*(1:96)/96)),2));
s=permute(s,[2 1 3]);

% imshow3(sos(s,3))
S0 = fif(mean(s(:,[1 3],:),2));
S1 = fif(mean(s(:,[2],:),2));
for cnt=1:size(S0,3)
    [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'}); 
end
k_gc = phzapply( k(:,:,:), y);
%imshow3(sos(k_gc,3))
% imshow3(sos(fft2c(k_gc),3))

