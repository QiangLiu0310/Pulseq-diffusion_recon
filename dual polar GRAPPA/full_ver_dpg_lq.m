%% full version of DPG from raw data Bruker
% Qiang Liu 10/Dec/2021
% MRLab SMU
clc; close all;clear all;
%% read and reconstruction k-space data
% shot 1 and shot 2
    path = ['D:\buda_data_21\20211228_190552_LQ_1_105\47']; 
    fid = fopen([path,'\fid'],'r');
    Rawdata = fread(fid,'int32');
    fclose(fid);
    CompData = Rawdata(1:2:end)+1i.*Rawdata(2:2:end);
    info = VT_BrukerInfo(path,1);
    % for multi-b value
    info.NDir=2+30*5;
    [kx,b0,PhaseCor] = kxtrajectory1(path,info.DIM(1));
    Prodata1 = reshape(CompData(1:end), [], info.Ncoil, info.DIM(3),info.NShot,info.NEX,info.NDir);%for data acquisition without navigator
    Prodata1 = reshape(Prodata1(end+1-info.DIM(1)*info.DIM(2)/info.NShot:end,:,:,:,:,:),info.DIM(1), info.DIM(2)/info.NShot, info.Ncoil, info.DIM(3),info.NShot,info.NEX,info.NDir);
    Prodata = permute(Prodata1, [2 1 3 4 5 6 7 8]);
    Prodata(2:2:end,:,:,:,:,:,:)=flipdim(Prodata(2:2:end,:,:,:,:,:),2);% flip the even readout lines since it is EPI.
    Prodata(:,:,:,:,2,:,:)=flipdim(Prodata(:,:,:,:,2,:,:),1);% flip the second shot for blip down acqusition.
    shot1=reshape(Prodata(:,:,:,:,1,:,:),[48 128 4 1 152]); % blip up acq
    shot2=reshape(Prodata(:,:,:,:,2,:,:),[48 128 4 1 152]); % blip down acq
% A0 data
    path = ['D:\buda_data_21\20211228_190552_LQ_1_105\43']; %D:\IRIS_Kidney\20210805_182105_LQ_1_99\14
    fid = fopen([path,'\fid'],'r');
    Rawdata = fread(fid,'int32');
    fclose(fid);
    CompData = Rawdata(1:2:end)+1i.*Rawdata(2:2:end);
    info = VT_BrukerInfo(path,1);
    [kx,b0,PhaseCor] = kxtrajectory1(path,info.DIM(1));
    Prodata1 = reshape(CompData(1:end), [], info.Ncoil, info.DIM(3),info.NShot,info.NDir);%for data acquisition without navigator
    Prodata1 = reshape(Prodata1(end+1-info.DIM(1)*info.DIM(2)/info.NShot:end,:,:,:,:,:),info.DIM(1), info.DIM(2)/info.NShot, info.Ncoil, info.DIM(3),info.NShot,info.NDir);
    pro1_2 = permute(Prodata1, [2 1 3 4 5 6 7 8]); 
% A1 data
    path = ['D:\buda_data_21\20211228_190552_LQ_1_105\44']; %D:\IRIS_Kidney\20210805_182105_LQ_1_99\16
    fid = fopen([path,'\fid'],'r');
    Rawdata = fread(fid,'int32');
    fclose(fid);
    CompData = Rawdata(1:2:end)+1i.*Rawdata(2:2:end);
    info = VT_BrukerInfo(path,1);
    [kx,b0,PhaseCor] = kxtrajectory1(path,info.DIM(1));
    Prodata1 = reshape(CompData(1:end), [], info.Ncoil, info.DIM(3),info.NShot,info.NDir);%for data acquisition without navigator
    Prodata1 = reshape(Prodata1(end+1-info.DIM(1)*info.DIM(2)/info.NShot:end,:,:,:,:,:),info.DIM(1), info.DIM(2)/info.NShot, info.Ncoil, info.DIM(3),info.NShot,info.NDir);
    pro3_4 = permute(Prodata1, [2 1 3 4 5 6 7 8]); 
% A2 data (negative)
    path = ['D:\buda_data_21\20211228_190552_LQ_1_105\45']; 
    fid = fopen([path,'\fid'],'r');
    Rawdata = fread(fid,'int32');
    fclose(fid);
    CompData = Rawdata(1:2:end)+1i.*Rawdata(2:2:end);
    info = VT_BrukerInfo(path,1);
    [kx,b0,PhaseCor] = kxtrajectory1(path,info.DIM(1));
    Prodata1 = reshape(CompData(1:end), [], info.Ncoil, info.DIM(3),info.NShot,info.NDir);%for data acquisition without navigator
    Prodata1 = reshape(Prodata1(end+1-info.DIM(1)*info.DIM(2)/info.NShot:end,:,:,:,:,:),info.DIM(1), info.DIM(2)/info.NShot, info.Ncoil, info.DIM(3),info.NShot,info.NDir);
    pro5_6 = permute(Prodata1, [2 1 3 4 5 6 7 8]); %  48    96     4     6     2     2
% A3 data
    path = ['D:\buda_data_21\20211228_190552_LQ_1_105\46']; 
    fid = fopen([path,'\fid'],'r');
    Rawdata = fread(fid,'int32');
    fclose(fid);
    CompData = Rawdata(1:2:end)+1i.*Rawdata(2:2:end);
    info = VT_BrukerInfo(path,1);
    [kx,b0,PhaseCor] = kxtrajectory1(path,info.DIM(1));
    Prodata1 = reshape(CompData(1:end), [], info.Ncoil, info.DIM(3),info.NShot,info.NDir);%for data acquisition without navigator
    Prodata1 = reshape(Prodata1(end+1-info.DIM(1)*info.DIM(2)/info.NShot:end,:,:,:,:,:),info.DIM(1), info.DIM(2)/info.NShot, info.Ncoil, info.DIM(3),info.NShot,info.NDir);
    pro7_8 = permute(Prodata1, [2 1 3 4 5 6 7 8]); %  48    96     4     6     2     2
    
%% dpg for blip down (readout+)
k=shot1(:,:,:,1,:);
% k=shot1(:,:,:,2,:);
% pro1_2=pro1_2(:,:,:,2,:);
% pro3_4=pro3_4(:,:,:,2,:);
% shot1=shot1(:,:,:,2,:);
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
%  imshow3(sos(fft2c(C),3))
%  imshow3(flipdim(permute(sos(fft2c(C),3),[2 1 3]),1))
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
Fdpg=zeros(96,128,4,size(shot1,5));
for i=1:size(shot1,5)
Fin{1} = zeros([96 128 4]); 
Fin{1}(2:4:end,:,:) =squeeze(k(1:2:end,:,:,:,i));
Fin{2} = zeros([96 128 4]);
Fin{2}(4:4:end,:,:) =squeeze(k(2:2:end,:,:,:,i));
Fdpg(:,:,:,i) = dpg_recon( Fin, N, 2, '5x9'); 
% Fdpgg(:,:,:,i)=Fdpg;
end% lq 2-'59'
% imshow3(sos(Fdpg),3))
imshow3(flipdim(permute(sos(fft2c(Fdpg),3),[2 1 3 4]),1))
blip_down=Fdpg;
%% dpg for blip up (readout-)
 k=shot2(:,:,:,1,:);
 k=flipdim(k,1);
% k=ifi(tmult(fif(k),diag(exp(-j*2*pi*(1:48)/48)),1));
%  a=k;
% k(1:2:end,:,:)=k(2:2:end,:,:);
% k(2:2:end,:,:)=a(1:2:end,:,:);
%  k(1:2:end,:,:)=ifi(tmult(fif(k(1:2:end,:,:)),diag(exp(-j*10.5*pi*(1:128)/128)),2));
%  k(2:2:end,:,:)=ifi(tmult(fif(k(2:2:end,:,:)),diag(exp(j*9.5*pi*(1:128)/128)),2));
% shift from 43 to 47
% you can test from code: imshow3(sos(pro1_2(1:2:end,:,:,:,1),3) sos(pro3_4(2:2:end,:,:,:,1),3))
A = zeros([96 128 4 1 ]);
A(1:4:end,:,:,:) = pro7_8(1:2:end,:,:,:,1);
A(2:4:end,:,:,:) = pro7_8(1:2:end,:,:,:,2);
A(3:4:end,:,:,:) = ifi(tmult(fif(pro5_6(2:2:end,:,:,:,1)),diag(exp(-j*7.6*pi*(1:128)/128)),2)); %7.6 is the best
A(4:4:end,:,:,:) = ifi(tmult(fif(pro5_6(2:2:end,:,:,:,2)),diag(exp(-j*7.6*pi*(1:128)/128)),2));
B = zeros([96 128 4 1 ]);
B(1:4:end,:,:,:) = ifi(tmult(fif(flipdim(pro5_6(1:2:end,:,:,:,1),2)),diag(exp(-j*7.6*pi*(1:128)/128)),2));
B(2:4:end,:,:,:) = ifi(tmult(fif(flipdim(pro5_6(1:2:end,:,:,:,2),2)),diag(exp(-j*7.6*pi*(1:128)/128)),2));
B(3:4:end,:,:,:) = flipdim(pro7_8(2:2:end,:,:,:,1),2);
B(4:4:end,:,:,:) = flipdim(pro7_8(2:2:end,:,:,:,2),2);
C = ifi(phzshift(  fif(A), fif(B), {'nofft'}));
kin{1} = A; 
kin{2} = B; 
kin{3} = C;
[~,~,N]=dpg_cal( size(kin{3}), kin, 1:size(kin{3},1), 'kernel','5x9','dks',2);
Fdpg=zeros(96,128,4,size(shot2,5));
for i=1:size(shot2,5)
Fin{1} = zeros([96 128 4]); 
Fin{1}(2:4:end,:,:) =squeeze(k(1:2:end,:,:,:,i));
Fin{2} = zeros([96 128 4]);
Fin{2}(4:4:end,:,:) =squeeze(k(2:2:end,:,:,:,i));
Fdpg(:,:,:,i) = dpg_recon( Fin, N, 2, '5x9'); 
end
 imshow3(flipdim(permute(sos(fft2c(Fdpg),3),[2 1 3 4]),1))
blip_up=Fdpg;

%% for FSL
blipd=squeeze(sos(fft2c(blip_down(:,:,:,1)),3));
blipu=squeeze(sos(fft2c(blip_up(:,:,:,1)),3));
blipd=flipdim(permute(repmat(blipd,1,1,1,4),[2 1 4 3]),1);
blipu=flipdim(flipdim(permute(repmat(blipu,1,1,1,4),[2 1 4 3]),1),2);

blipdnii= make_nii(blipd(:,:,:,1));
save_nii(blipdnii,'blip_down.nii')
blipunii= make_nii(blipu(:,:,:,1));
save_nii(blipunii,'blip_up.nii')