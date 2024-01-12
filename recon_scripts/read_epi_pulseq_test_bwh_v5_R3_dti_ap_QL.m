% Pulseq-Siemens data reconstruction
% using mapVBVD https://github.com/CIC-methods/FID-A/tree/master/inputOutput/mapVBVD to read the data
% Pulseq toolbox is required to calculate the trajectory for EPI
%% EPI recon includes: 
% 1. ramp sampling intepolaration (Pulseq code) 
% 2. Ghost correction (LPC based on Dr.William Scott Hoge's code) 
% 3. GRAPPA reconstruction
% 4. Homodyne for Partial Fourier recon
% 5. Apodization (based on Dr.Bilgic Berkin's toolbox @Pulseq_Mprage_Recon_Toolbox)
%% 
% by Qiang Liu et al., Brigham and Women's Hospital, Harvard Medical School
% qliu30@mgh.harvard.edu
% qiangliu2019@163.com

%%
close all; clear all ; clc
% addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Pulseq_Mprage_Recon_Toolbox'))
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Pulseq-diffusion'))

%%
data_path = '/data/pnl/home/ql087/data_bwh/2023_10_04_bwh_pulseq_std_sub1/'
filename1 = 'meas_MID00673_FID959057_pulseq_std_88sli_30dir_B'; % image data
filename2 = 'meas_MID00672_FID959056_pulseq_std_88sli_ref'; % reference data

D=dir([data_path filename1]);
[~,I]=sort([D(:).datenum]);
twix_obj = mapVBVD([data_path filename1]);
seq = mr.Sequence();
read(seq,'/Users/ln915/Documents/code/github_code/Pulseq_recon_scripts/seq_file/epidiff_3_shot_2us_1p5_88sli_30dir_1_appa_recon.seq')
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

%% define raw data

if iscell(twix_obj)
    rawdata = double(twix_obj{end}.image.unsorted());
else
    rawdata = double(twix_obj.image.unsorted());
end
rawdata=single(rawdata);
traj_recon_delay=0e-6;
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
%% automatic detection of the measurement parameters (FOV, matrix size, etc)
nADC = size(rawdata, 1);
k_last=ktraj_adc(:,end);
k_2last=ktraj_adc(:,end-nADC);
delta_ky=k_last(2)-k_2last(2);
fov=1/abs(delta_ky);

%% manually (re-)define FOV and resolution
Nx=146; Ny=Nx;
R=3;
Rseg=1;
pf=6/8; 
fe_pf=1;
Nx=ceil(Nx*fe_pf); 
Ny_sampled=ceil(Ny*pf/R/Rseg);
slice_num=88;
%% analyze the trajecotory, resample the data
nCoils = size(rawdata, 2); % the incoming data order is [kx coils acquisitions]
nAcq=size(rawdata,3);
nD=size(ktraj_adc, 1);

kxmin=min(ktraj_adc(1,:));
kxmax=max(ktraj_adc(1,:));
kxmax1=kxmax/(Nx/2-1)*(Nx/2); % this compensates for the non-symmetric center definition in FFT
kmaxabs=max(kxmax1, -kxmin);

kxx= ((-Nx/2):(Nx/2-1))/(Nx/2)*kmaxabs; % kx-sample positions

ktraj_adc2=reshape(ktraj_adc,[size(ktraj_adc,1), nADC, size(ktraj_adc,2)/nADC]);
t_adc2=reshape(t_adc,[nADC, length(t_adc)/nADC]);

diff_plus_dum=33; % b0 (dummy) + b0 +b0 (PA)+ b1000
nAd=nAcq/diff_plus_dum;

data_resampled=zeros(length(kxx),nCoils, nAd, diff_plus_dum);
ktraj_adc2=reshape(ktraj_adc2,[size(ktraj_adc2,1) size(ktraj_adc2,2) nAd 3]); 
rawdata=reshape(rawdata,[size(rawdata,1)  size(rawdata,2) nAd diff_plus_dum]);

delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object

tic
parfor t_c=1:size(data_resampled,4)
    for a=1:nAd
        for c=1:nCoils
            data_resampled(:,c,a,t_c)=interp1(ktraj_adc2(1,:,a),rawdata(:,c,a, t_c),kxx,'spline',0);
        end
    end
end
toc
delete(gcp('nocreate'))
disp('interporlation completed for std data')
data_resampled=reshape(data_resampled,[length(kxx), nCoils, nAcq]);

%% save the reference data
ref=permute(data_resampled,[1 3 2]); 
ref=reshape(ref,[Nx Ny_sampled slice_num*Rseg*diff_plus_dum nCoils]);
multiply_ref=slice_num*Rseg;

%% start to do LPC for EPI
k=ref(:,:,(multiply_ref+1):end,:); % k-space data
k=reshape(k,[Nx Ny_sampled multiply_ref diff_plus_dum-1 nCoils]);
k_gc=zeros(Ny_sampled,Nx,multiply_ref, (diff_plus_dum-1),nCoils);
ref=ref(:,:,1:multiply_ref,:);

for seg_loop=1:(diff_plus_dum-1)

    for i=1:size(ref,3)

        s=squeeze(ref(:,2:4,i,:)); %pick the center three lines will be better
        k_1=squeeze(k(:,:,i,seg_loop,:));
        % LPC based on WSH's LPC code
        S0 = fif(mean(s(:,[1 3],:),2));
        S1 = fif(mean(s(:,[ 2 ],:),2));
        for cnt=1:size(S0,3); [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'}); end;
        k_gc(:,:,i,seg_loop,:) = phzapply( permute(k_1(:,1:end,:),[2 1 3]), y);

    end
end


clear rawdata data_resampled
clear t_adc t_adc2 t_ktraj ktraj ktraj_adc ktraj_adc2 s S0 S1 data_resampled

k_gc=squeeze(reshape(k_gc,[Ny_sampled, Nx, slice_num,Rseg , (diff_plus_dum-1), nCoils]));
k_gc=permute(k_gc,[2 5 1 3 4]);
sz_k=[Nx nCoils Ny_sampled slice_num diff_plus_dum-1];
Kimage_first = reshape(k_gc,sz_k);
Kimage_first_full=zeros(Nx, nCoils, Ny_sampled*3,slice_num,diff_plus_dum-1);
Kimage_first_full(:,:,R:R:end,:,:) = Kimage_first;
clear rawdata data_resampled  Kimage_first k_gc

%% ACS data
D=dir([data_path filename2]);
[~,I]=sort([D(:).datenum]);
twix_obj = mapVBVD([data_path filename2]);
seq = mr.Sequence();
read(seq,'/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Tests/Test_49_May_23/recon/epidiff_3_shot_ref_1p5mm_88sli_appa_recon.seq') % I want to keep this line, QL
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

%% define raw data

if iscell(twix_obj)
    rawdata = double(twix_obj{end}.image.unsorted());
else
    rawdata = double(twix_obj.image.unsorted());
end
rawdata=single(rawdata);
rawdata=rawdata(:,:,1:size(rawdata,3)*2/3); %extract the AP data only
traj_recon_delay=0e-6;
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);
%% automatic detection of the measurement parameters (FOV, matrix size, etc)
nADC = size(rawdata, 1);
k_last=ktraj_adc(:,end);
k_2last=ktraj_adc(:,end-nADC);
delta_ky=k_last(2)-k_2last(2);
fov=1/abs(delta_ky);

%% manually (re-)define FOV and resolution
Nx=146; Ny=Nx; Ny_sampled=ceil(Ny*6/8/3);
%% analyze the trajecotory, resample the data
nCoils = size(rawdata, 2); % the incoming data order is [kx coils acquisitions]
nAcq=size(rawdata,3);
nD=size(ktraj_adc, 1);

kxmin=min(ktraj_adc(1,:));
kxmax=max(ktraj_adc(1,:));
kxmax1=kxmax/(Nx/2-1)*(Nx/2); % this compensates for the non-symmetric center definition in FFT
kmaxabs=max(kxmax1, -kxmin);

kxx= ((-Nx/2):(Nx/2-1))/(Nx/2)*kmaxabs; % kx-sample positions
ktraj_adc2=reshape(ktraj_adc,[size(ktraj_adc,1), nADC, size(ktraj_adc,2)/nADC]);
t_adc2=reshape(t_adc,[nADC, length(t_adc)/nADC]);

data_resampled=zeros(length(kxx), nCoils, nAcq);
for a=1:nAcq
    for c=1:nCoils
        data_resampled(:,c,a)=interp1(ktraj_adc2(1,:,a),rawdata(:,c,a),kxx,'spline',0);
    end
end

ref=permute(data_resampled,[1 3 2]);
diff_plus_dum=2;
ref=reshape(ref,[Nx Ny_sampled slice_num*R*(diff_plus_dum) nCoils]);
multiply_ref=slice_num*R;

%% start to do LPC for EPI
k=ref(:,:,(multiply_ref+1):end,:);
k_gc_1=zeros(Ny_sampled,Nx,multiply_ref,nCoils);
ref=ref(:,:,1:multiply_ref,:);

for i=1:size(ref,3)

    s=squeeze(ref(:,2:4,i,:));
    k_1=squeeze(k(:,:,i,:));
    %% LPC based on WSH's LPC code
    S0 = fif(mean(s(:,[1 3],:),2));
    S1 = fif(mean(s(:,[ 2 ],:),2));
    for cnt=1:size(S0,3); [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'}); end;
    k_gc_1(:,:,i,:) = phzapply( permute(k_1(:,1:end,:),[2 1 3]), y);

end

clear rawdata data_resampled

acs=permute(k_gc_1(7:18,:,:,:), [2 4 1 3]); % use the central 12 lines (please check the phase of the image to make sure it is in the center)
acs=reshape(acs,[Nx nCoils 12 slice_num R]);
refscan_first=zeros(Nx,nCoils,12*R,slice_num);
for iii=1:R
    refscan_first(:,:,iii:R:end,:)=squeeze(acs(:,:,:,:,iii));
end

clear k_gc_1 t_adc t_adc2 t_ktraj ktraj ktraj_adc ktraj_adc2 rawdata S0 S1 s ref

% GRAPPA Kernel

if R~=1

    if ~exist('chi','var'); chi = 1e-6; end;
    if ~exist('eta','var'); eta = 1; end;
    %  the next step is to generate the grappa parameters

    if R~=1
        for slc = 1:slice_num
            [~,~,Ng_first{slc}]=recongrappa_multik([36 Nx nCoils],permute(refscan_first(:,:,:,slc),[3 1 2]),[],'kernel','2x5','dks',R*[1 2],...
                'chi',chi,'eta',eta );
        end
    end

end

% change the loop structure
Kimage_first_full=reshape(Kimage_first_full,[size(Kimage_first_full,1), size(Kimage_first_full,2), size(Kimage_first_full,3), size(Kimage_first_full,4)*size(Kimage_first_full,5)]);
Ng_first_1=repmat(Ng_first,[1 1  size(Kimage_first_full,4)/size(Ng_first,2)]);
Ng_first_1=reshape(Ng_first_1,[size(Ng_first_1,1), size(Kimage_first_full,4)]);
I_short=zeros(size(Kimage_first_full,1),size(Kimage_first_full,1),size(Kimage_first_full,4));
img_coil1=zeros(Nx,Ny,nCoils);

for slc = 1:size(Kimage_first_full,4)
    img_temp1= recongrappa_multik([Ny_sampled*R Nx nCoils],permute(Kimage_first_full(:,:,:,slc),[3 1 2]),[],'kernel','2x5','N',Ng_first_1{slc});
    % correct the partial-Fourier acquiisition, for each coil
    for cntC=1:nCoils
        img_coil1(:,:,cntC) = flip(reconhd( img_temp1(:,:,cntC),size(img_temp1,1),size(img_temp1,2) )) ;
        img_coil1(:, :,cntC)= apodize3d_QL_v2(img_coil1(:,:,cntC), 0.25);
    end
    I_short(:,:,slc) = sqrt(sum(abs( fif( img_coil1 ) ).^2,3));
    slc
end

I_short=reshape(I_short, [Nx Ny slice_num 32]);
p=[2:2:slice_num, 1:2:slice_num];
[~,or]=sort(p);
I_short=I_short(:,:,or,:);
I_short=I_short(:,:,end:-1:1,:);

% save_path='/data/pnl/home/ql087/data_processing/2023_12_08_pulseq_multiecho/';
% save([save_path 'multiecho_te1.mat'], 'I_short', '-v7.3')

