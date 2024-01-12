%--------------------------------------------------------------------------
%% load mprage: 
%--------------------------------------------------------------------------

addpath utils/

filename = 'meas_MID00113_FID82601_pulseq_140.dat';

% dimensions:
% 1024 x 36 x 27840 = 1280 x 36 x (192x(144+1))
% 4x ro oversampled
% 36 chan, 144 PE lines + 1 for noise reference

tic
    dt = mapVBVD(filename);
    
    res = dt{end}.image.unsorted();     
toc



%--------------------------------------------------------------------------
%% define GRAPPA undersampling pattern along slowest dimension (outer pe loop)
%--------------------------------------------------------------------------

Nx = 256;
Ny = 256;
Nz = 192;

    
N = [192 256 256];                % matrix sizes

ax=struct; % encoding axes
ax.d1='z'; % the fastest dimension (readout)
ax.d2='x'; % the second-fast dimension (the inner pe loop)
ax.d3=setdiff('xyz',[ax.d1 ax.d2]); % automatically set the slowest dimension
ax.n1=strfind('xyz',ax.d1);
ax.n2=strfind('xyz',ax.d2);
ax.n3=strfind('xyz',ax.d3);

nACS = 32;  % number of auto-calibration lines (dense sampling in center)
R = 2;      % undersampling outside ACS region
S = 0*(1:N(ax.n3));
S(1:R:end) = 1;
S( (end/2-nACS/2):(end/2+nACS/2-1) ) = 1;  % TODO: uncomment
J = find(S == 1);

num_chan = s(res,2);
num_ro = s(res,1);

mask = zeross([Nz, Ny]);
mask(:,J) = 1;

%--------------------------------------------------------------------------
% sampling mask
%--------------------------------------------------------------------------

m2d = zeros([Nz, Ny]);

m2d( :, 1:2:end ) = 1;
m2d(:, end/2-nACS/2 : end/2+nACS/2-1) = 1;

mosaic(m2d,1,1,1)

sum(sum(m2d,1)~=0)  % number of sampled lines

    

%--------------------------------------------------------------------------
%% crop readout
%--------------------------------------------------------------------------

kspace = ifftc(res,1);
kspace = fftc(kspace(1+end/2-Nx/2:end/2+Nx/2,:,:,:),1);


%--------------------------------------------------------------------------
% last TR is for noise scans: 
%--------------------------------------------------------------------------

kspace_temp = kspace(:,:,1:end-Nz);
kspace_temp = reshape(kspace_temp, [Nx, num_chan, Nz, length(J)]);

kspace_pad = zeross([Nx, num_chan, Nz, Ny]);
kspace_pad(:,:,:,J) = reshape(kspace_temp, [Nx, num_chan, Nz, length(J)]);
kspace_pad = permute(kspace_pad, [1,3,4,2]);    


% reference data
kspace_patref = zeross(size(kspace_pad));
kspace_patref(:,:,end/2-nACS/2:end/2+nACS/2-1,:) = kspace_pad(:,:,end/2-nACS/2:end/2+nACS/2-1,:);


% go to hybrid space: ifft over kx
img_hybrid = ifftc(kspace_pad,1);               %x,ky,kz,chan
img_hybrid_patref = ifftc(kspace_patref,1);               %x,ky,kz,chan


clear res


%--------------------------------------------------------------------------
%% scale data: make patref image unit norm
%--------------------------------------------------------------------------

scale_factor = 1 / norm2(img_hybrid_patref)

img_hybrid = img_hybrid * scale_factor;
img_hybrid_patref = img_hybrid_patref * scale_factor;

imagesc3d2(rsos(ifftc(ifftc(img_hybrid_patref,3),2),4), s(img_hybrid_patref)/2, 1, [90,180,180], [0.,7.5e-4]),



%--------------------------------------------------------------------------
%% esprit coil sensitivity estimation: 
%% can skip tp next cell and load pre-computed maps
%--------------------------------------------------------------------------

num_acs = 24;
kernel_size = [6,6];
eigen_thresh = 0.5;     % get large coil sens mask since we just need to combine coils and not recon

delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object

total_cores = c.NumWorkers; % use half the cores
parpool(ceil(total_cores/2))

receive = zeross(size(img_hybrid_patref));

tic    
parfor slc_select = 1:Nx
    disp(num2str(slc_select))
    
    kspace_slice = sq(img_hybrid(slc_select,:,:,:));
    kspace_slice = circshift(kspace_slice, [0,1]);

    Kspace_Acs = zeross(size(kspace_slice));
    Kspace_Acs(:,1+end/2-nACS/2:end/2+nACS/2,:) = kspace_slice(:,1+end/2-nACS/2:end/2+nACS/2,:);
    
    [maps, weights] = ecalib_soft( Kspace_Acs, num_acs, kernel_size, eigen_thresh );

    receive(slc_select,:,:,:) = dot_mult(maps, weights >= eigen_thresh); 
end
toc

delete(gcp('nocreate'))


%--------------------------------------------------------------------------
%% load pre-computed coil sensitivities
%--------------------------------------------------------------------------

load receive


%--------------------------------------------------------------------------
%% loop over slices
%--------------------------------------------------------------------------
 
delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object

total_cores = c.NumWorkers; % 48 cores for marduk
parpool(ceil(total_cores/2))

num_eco = 1;
substitute_acs = 0;

Rx = 1;
Ry = 2;

delx = zeros(num_chan,1);        % starting index of ky lines
dely = ones(num_chan,1);

lambda_percent = 10;             % percentage of sigma_min to use as tikhonov regularizer

Img_Res = zeross([Nx,Nz,Ny,num_chan]);
eN = [Nz,Ny];


tic
parfor slc_select = 1:Nx
    disp(num2str(slc_select))
    
    
    kspace_slice = sq(img_hybrid(slc_select,:,:,:));
    kspace_slice = circshift(kspace_slice, [0,1]);
    
    img_slice = ifft2call(kspace_slice);

    Kspace_Sampled = kspace_slice;

    Kspace_Acs = zeross(size(Kspace_Sampled));
    Kspace_Acs(:,1+end/2-nACS/2:end/2+nACS/2,:) = kspace_slice(:,1+end/2-nACS/2:end/2+nACS/2,:);
    
    
    num_acs = [eN(1), nACS]-2;        % size reduced due to 1 voxel circshift
    kernel_size = [3,3];              % odd kernel size

    Img_Grappa = zeros([eN(1:2), num_chan, num_eco]);

    for t = 1:num_eco
        Img_Grappa(:,:,:,t) = grappa_gfactor_2d_jvc3( Kspace_Sampled(:,:,:,t), Kspace_Acs, Rx, Ry, num_acs, kernel_size, lambda_percent, substitute_acs, delx, dely );
    end
    
    temp = fft2call(Img_Grappa);
    temp(:,1+end/2-nACS/2:end/2+nACS/2,:) = Kspace_Acs(:,1+end/2-nACS/2:end/2+nACS/2,:);
    Img_Grappa = ifft2call(temp);    
    
    Img_Res(slc_select,:,:,:) = Img_Grappa; 
end
toc

delete(gcp('nocreate'))


imagesc3d2(rsos(Img_Res,4), s(Img_Res)/2, 22, [-90,180,180], [0.,7.5e-4]),


%--------------------------------------------------------------------------
%% apodize, coil combine 
%--------------------------------------------------------------------------

Img_Res_apo0p25 = apodize3d(Img_Res, 0.25);

Img_espirit_apo0p25 = coil_combine(Img_Res_apo0p25, receive, 4);


imagesc3d2(Img_espirit_apo0p25, s(Img_espirit_apo0p25)/2, 23, [-90,180,180], [0.,7.5e-4]),

 