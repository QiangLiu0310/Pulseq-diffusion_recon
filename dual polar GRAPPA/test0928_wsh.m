%% DPG for BUDA 
%% by Scott Hoge
%% Qiang Liu
%% 28th/Sep/2021
%% MR Lab
%% for blip up

clc;clear all; close all;
% load the calibration data
load('kidneypro1_2.mat')
load('kidneypro3_4.mat')
%  load('pro1_2.mat')
%  load('pro3_4.mat')
load('kidneyshot1.mat')
%  load('shot1.mat')
%  load('noblip.mat')
 k=shot1(:,:,:,1,1);


% interleave to form the RO+ image
% shift from 43 to 47
% you can test from code: imshow3(sos(pro1_2(1:2:end,:,:,:,1),3) sos(pro3_4(2:2:end,:,:,:,1),3))
A = zeros([96 128 4 1 ]);
A(1:4:end,:,:,:) = pro1_2(1:2:end,:,:,:,1);
A(2:4:end,:,:,:) = pro1_2(1:2:end,:,:,:,2);
% A(3:4:end,:,:,:) = pro3_4(2:2:end,:,:,:,1);
% A(4:4:end,:,:,:) = pro3_4(2:2:end,:,:,:,2);
%  imshow3([sos(pro1_2(1:2:end,:,:,:,1),3); sos(pro3_4(2:2:end,:,:,:,1),3)])
A(3:4:end,:,:,:) = ifi(tmult(fif(pro3_4(2:2:end,:,:,:,1)),diag(exp(-j*7.6*pi*(1:128)/128)),2)); %7.6 is the best
A(4:4:end,:,:,:) = ifi(tmult(fif(pro3_4(2:2:end,:,:,:,2)),diag(exp(-j*7.6*pi*(1:128)/128)),2));
% imshow3(sos(fft2c(A),3))
% imshow3(flipdim(permute(sos(fft2c(A),3),[2 1 3]),1))
% A=A(:,:,:,1);
% % interleave to form the RO- image
B = zeros([96 128 4 1 ]);
% B(1:4:end,:,:,:) = pro3_4(1:2:end,:,:,:,1);
% B(2:4:end,:,:,:) = pro3_4(1:2:end,:,:,:,2);
B(1:4:end,:,:,:) = ifi(tmult(fif(flipdim(pro3_4(1:2:end,:,:,:,1),2)),diag(exp(-j*7.6*pi*(1:128)/128)),2));
B(2:4:end,:,:,:) = ifi(tmult(fif(flipdim(pro3_4(1:2:end,:,:,:,2),2)),diag(exp(-j*7.6*pi*(1:128)/128)),2));
B(3:4:end,:,:,:) = flipdim(pro1_2(2:2:end,:,:,:,1),2);
B(4:4:end,:,:,:) = flipdim(pro1_2(2:2:end,:,:,:,2),2);
% imshow3(sos(fft2c(B),3))
% B=B(:,:,:,1);
% imshow3(flipdim(permute(sos(fft2c(B),3),[2 1 3]),1))


C = ifi(phzshift(  fif(A), fif(B), {'nofft'}));
 imshow3(sos(fft2c(C),3))
 imshow3(flipdim(permute(sos(fft2c(C),3),[2 1 3]),1))
% % B = zeros(size(acs1));
% % B(:,1:2:end,:) = acs2(:,1:2:end,:);
% % B(:,2:2:end,:) = acs1(:,2:2:end,:);
% % B = permute(B,[2 1 3]);

% combine coherently, to form the PLACE image
% the lines below calls a MEX function from the Fast Imaging Library 
% C = ifi(pos_neg_add( fif(A),fif(B) ))/2;


% if the FIL function is not available, you can replace the above line(s) with
% the call below, to read in the result from a file.
% C = readnd('dpg_C.nd');

% calibrate the DPG coefficients
kin{1} = A; 
kin{2} = B; 
kin{3} = C;
[~,~,N]=dpg_cal( size(kin{3}), kin, 1:size(kin{3},1), 'kernel','5x9','dks',2);

% read in the raw EPI data (R=2)
% k = readnd('dpg_k.nd');
% reconstruct the data using DPG
% k=ifi(tmult(fif(k),diag(exp(-j*1*pi*(1:96)/96)),2));

Fin{1} = zeros([96 128 4]); 
Fin{1}(2:4:end,:,:) =k(1:2:end,:,:);
Fin{2} = zeros([96 128 4]);
Fin{2}(4:4:end,:,:) =k(2:2:end,:,:);
Fdpg = dpg_recon( Fin, N, 2, '5x9'); % lq 2-'59'
% imshow3(sos(fft2c(Fdpg),3))
imshow3(flipdim(permute(sos(fft2c(Fdpg),3),[2 1 3]),1))
% compare with an LPC+GRAPPA recon:
c=zeros(48,96,4);
[~,~,Ng]=recongrappa(size(kin{3}),kin{3},vec(1:96),'kernel','5x9','dks',[2;4]); %96

% compute the standard Nyquist ghost correction coefficients:
% s = readnd('dpg_k_pcref.nd');
% S0 = fif(mean(s(:,[1 3],:),2));
% S1 = fif(mean(s(:,[ 2 ],:),2));
% for cnt=1:size(S0,3); [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'}); end;
% % y = readnd('y.nd');                     % derived from dpg_k_pcref.nd
% 
% k_gc = phzapply( permute(k(:,1:2:end,:),[2 1 3]), y);


% compute the standard Nyquist ghost correction coefficients:
s2 = noblip;
s2(2:2:end,:,:) = flipdim(s2(2:2:end,:,:),2);
s=s2(46:48,:,:,1,1);

% s= ifi(tmult(fif(s),diag(exp(j*4*pi*(1:96)/96)),2));
s=permute(s,[2 1 3]);

% imshow3(sos(s,3))
S0 = fif(mean(s(:,[1 3],:),2));
S1 = fif(mean(s(:,[2],:),2));
for cnt=1:size(S0,3)
    [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'}); 
end
k_gc = phzapply( k(:,:,:), y);

Flpc = recongrappa([96 96 4],k_gc,vec(1:2:96),'kernel','5x9','dks',2,'N',Ng);

imagesc(sqrt(sum(abs([ flipdim(fif(Flpc),1); flipdim(fif(Fdpg),1) ]).^2,3))); %32 channels, ql
axis('image');
colormap(gray(256))