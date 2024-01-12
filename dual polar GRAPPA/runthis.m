clc;clear all; close all;
% load the calibration data
acs1 = readnd('dpg_acs1.nd');               % RO+ lead data
acs2 = readnd('dpg_acs2.nd');               % RO- lead data

% interleave to form the RO+ image
A = zeros(size(acs1));
A(:,1:2:end,:) = acs1(:,1:2:end,:);
A(:,2:2:end,:) = acs2(:,2:2:end,:);
A = permute(A,[2 1 3]);

% interleave to form the RO- image
B = zeros(size(acs1));
B(:,1:2:end,:) = acs2(:,1:2:end,:);
B(:,2:2:end,:) = acs1(:,2:2:end,:);
B = permute(B,[2 1 3]);

% combine coherently, to form the PLACE image
% the lines below calls a MEX function from the Fast Imaging Library 
% C = ifi(pos_neg_add( fif(A),fif(B) ))/2;
C = ifi(phzshift(  fif(A), fif(B), {'nofft'}));

% if the FIL function is not available, you can replace the above line(s) with
% the call below, to read in the result from a file.
% C = readnd('dpg_C.nd');

% calibrate the DPG coefficients
kin{1} = A; 
kin{2} = B; 
kin{3} = C;
[~,~,N]=dpg_cal( size(kin{3}), kin, 1:size(kin{3},1), 'kernel','2x5','dks',2);

% read in the raw EPI data (R=2)
k = readnd('dpg_k.nd');
% reconstruct the data using DPG
Fin{1} = zeros([60 96 32]); 
Fin{1}(2:4:end,:,:) = permute(k(:,1:4:end,:),[2 1 3]);
Fin{2} = zeros([60 96 32]);
Fin{2}(4:4:end,:,:) = permute(k(:,3:4:end,:),[2 1 3]);
Fdpg = dpg_recon( Fin, N, 2, 2);

% compare with an LPC+GRAPPA recon:
[~,~,Ng]=recongrappa(size(kin{3}),kin{3},vec(1:30),'kernel','2x5','dks',[2;4]); % 30 96 32/

% compute the standard Nyquist ghost correction coefficients:
s = readnd('dpg_k_pcref.nd');
S0 = fif(mean(s(:,[1 3],:),2));
S1 = fif(mean(s(:,[ 2 ],:),2));
for cnt=1:size(S0,3); [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'}); end;
% y = readnd('y.nd');                     % derived from dpg_k_pcref.nd

k_gc = phzapply( permute(k(:,1:2:end,:),[2 1 3]), y);
Flpc = recongrappa([60 96 32],k_gc,vec(1:2:60),'kernel','2x5','dks',2,'N',Ng);

imagesc(sqrt(sum(abs([ flipdim(fif(Flpc),1); flipdim(fif(Fdpg),1) ]).^2,3))); %32 channels, lq
axis('image');
colormap(gray(256))
