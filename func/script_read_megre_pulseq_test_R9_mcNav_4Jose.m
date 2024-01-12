%--------------------------------------------------------------------------
%% load 3d-gre: R1, 4mm
%--------------------------------------------------------------------------

addpath func/

data_path = [pwd, '/data/'];

    
% filename1 = 'meas_MID00880_FID119945_pulseq_140_gre_3mm_R1.dat';          % mtx 85x64x75
% filename1 = 'meas_MID00882_FID119947_pulseq_140_gre_3mm_R9.dat';          % mtx 85x72x81
% filename1 = 'meas_MID00884_FID119949_pulseq_140_gre_3mm_R9_varEs.dat';    % mtx 85x72x81
filename1 = 'meas_MID00883_FID119948_pulseq_140_gre_3mm_R9_mcNav.dat';    % mtx 85x72x81

tic
    dt1 = mapVBVD([data_path, filename1]);
    res1 = dt1{end}.image.unsorted();
toc



%--------------------------------------------------------------------------
%% sort data and load sequence params
%--------------------------------------------------------------------------

save_path = [pwd, '/20220421/'];

% filename = 'gre_3d_3mm_R1';
% filename = 'gre_3d_3mm_R9';
% filename = 'gre_3d_3mm_R9_varES';
filename = 'gre_3d_3mm_R9_mcNav';


load([save_path, 'SeqParamsUpdated_', filename, '.mat'])


num_chan = size(res1,2);
num_echoes = length(SeqParamsUpdated.TE);
Dims = SeqParamsUpdated.Dims;

Data = zeros([Dims,num_chan,num_echoes]);
NavFully = zeros([Dims(1:3),num_chan,num_echoes]);


Ny = Dims(2);
Nz = Dims(3);
nTE = num_echoes;

KspaceMaskKyKzt = SeqParamsUpdated.KspaceMaskKyKzt;

if isfield(SeqParamsUpdated, 'KspaceOrderNavKyKzt')
    disp('using mcNav')
    KspaceOrderNavKyKzt = SeqParamsUpdated.KspaceOrderNavKyKzt;
end

iShotData = 0;
iShotNav = 0;
Acq = 0;

for iShot = 1:length(SeqParamsUpdated.labelData0_Nav1)
    if ~(mod(iShot,100))
        disp([num2str(iShot), ' / ', num2str(length(SeqParamsUpdated.labelData0_Nav1))])
    end
    
    if SeqParamsUpdated.labelData0_Nav1(iShot) == 0
        iShotData = iShotData + 1;
        [PE1, PE2, TE] = ind2sub([Ny,Nz,nTE],find(KspaceMaskKyKzt==iShotData));

        for t = 1:num_echoes
            Acq = Acq + 1;
            Data(:, PE1(t), PE2(t),:,t) = res1(:,:,Acq);
        end
        
    else
        iShotNav = iShotNav + 1;
        [PE1, PE2, TE] = ind2sub([Ny,Nz,nTE],find(KspaceOrderNavKyKzt==iShotNav));
        
        for t = 1:num_echoes
            Acq = Acq + 1;

            NavFully(:, PE1(t), PE2(t),:,t) = res1(:,:,Acq);
        end

    end
end



%--------------------------------------------------------------------------
%% parse data
%--------------------------------------------------------------------------


img = ifft3call(Data);

% remove 2x oversampling
img = img(1+end/4:3*end/4,:,:,:,:);

img_R9 = img;


if isfield(SeqParamsUpdated, 'KspaceOrderNavKyKzt')
    img_nav = ifft3call(NavFully);

    img_nav = img_nav(1+end/4:3*end/4,:,:,:,:);

    for t = 1:num_echoes
        imagesc3d2(rsos(img_nav(:,:,:,:,t),4), s(img_nav)/2, 100+t, [90,90,-90], [0.,5e-3]), setGcf(.5)
        imagesc3d2(rsos(NavFully(:,:,:,:,t),4).^.5, s(img_nav)/2, 120+t, [0,0,0], [0.,1e-1]), setGcf(.5)
    end
end


for t = 1:num_echoes
    imagesc3d2(rsos(img(:,:,:,:,t),4), s(img)/2, t, [90,90,-90], [0.,5e-3]), setGcf(.5)
    imagesc3d2(rsos(Data(:,:,:,:,t),4).^.5, s(img)/2, 20+t, [0,0,0], [0.,1e-1]), setGcf(.5)
end



%--------------------------------------------------------------------------
%% load coil maps
%--------------------------------------------------------------------------

load([pwd, '/data/receive_3mm_R9.mat'])


%--------------------------------------------------------------------------
%% lsqr sense recon
%--------------------------------------------------------------------------


lsqr_iter = 200;
lsqr_tol = 1e-3;

img_res = zeross(s(sq(img_R9(:,:,:,1,:))));
m2d = sq(Data(1,:,:,:,:)) ~=0;

mosaic(sq(m2d(:,:,1,:)),2,4,56,'',[0,1]),colormap parula

    
delete(gcp('nocreate'))
c = parcluster('local');    % build the 'local' cluster object

total_cores = c.NumWorkers; % 48 cores for marduk
parpool(ceil(total_cores/2))

tic
parfor slc_select = 1:s(img_R9,1)
    disp(num2str(slc_select))
    
    sens = sq(receive(slc_select,:,:,:));
    kspace_slc = fft2call(sq(img_R9(slc_select,:,:,:,:)));


    param = [];
    param.N = Dims(2:3);
    param.num_chan = num_chan;
    param.lambda = 1e-4;        % L2 reg

    param.sens = sens;

    Res = zeross([param.N, num_echoes]);

    for t = 1:num_echoes
        kspace_coils = kspace_slc(:,:,:,t);
        param.m2d = m2d(:,:,:,t);

        res = lsqr(@apply_sense_tikc, cat(1, kspace_coils(:), zeros(prod(param.N),1)), lsqr_tol, lsqr_iter, [], [], [], param);  
        
        Res(:,:,t) = reshape(res, param.N);
    end

    img_res(slc_select,:,:,:) = Res;
end
toc

delete(gcp('nocreate'))
 


for t = 1:num_echoes
    imagesc3d2(img_res(:,:,:,t), s(img_res)/2, t, [90,90,-90], [0.,4e-3]), setGcf(.5)
end


% coil combine nav echoes
img_nav_combo = zeross(size(img_res));

for t = 1:num_echoes
    img_nav_combo(:,:,:,t) = coil_combine(img_nav(:,:,:,:,t), receive);
    
    imagesc3d2(img_nav_combo(:,:,:,t), s(img_res)/2, t+100, [90,90,-90], [0.,4e-3]), setGcf(.5)
end






