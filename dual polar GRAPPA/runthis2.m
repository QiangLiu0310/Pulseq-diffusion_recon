% an example of DPG w/o PLACE calibration data

acs2 = readnd('dpg_acs2.nd');           % RO- lead data
y2 = readnd('y_acs2.nd');               % ghost correction coefficients

% ghost-corrected ACS data
c2_gc = permute( phzapply( permute(acs2,[2 1 3]), y2),[2 1 3]);

% generate the GRAPPA coefficients 
[~,~,Ng2]=grappa(size(c2_gc),c2_gc,vec(1:30),'2x5',2);

% synthesize the RO+ and RO- source data from the ACS data using GRAPPA
acs2n = grappa(size(c2_gc),acs2(:,1:2:end,:),vec(1:2:30),'2x5',2,Ng2);
acs2p = grappa(size(c2_gc),acs2(:,2:2:end,:),vec(2:2:30),'2x5',2,Ng2);
acs2c = ifi(pos_neg_add( fif(acs2n),fif(acs2p) ))/2;

% calibrate DPG
kin{1} = permute(acs2p,[2 1 3]);        % RO+ in the first cell
kin{2} = permute(acs2n,[2 1 3]);        % RO- in the second cell
kin{3} = permute(acs2c,[2 1 3]);
[~,~,N2]=dpg_cal( size(kin{3}), kin, 1:size(kin{3},1), 'kernel','2x5','dks',2);


% reconstruct the accelerated data using DPG
k = readnd('dpg_k.nd');

Fin{1} = zeros([60 96 32]); Fin{1}(2:4:end,:,:) = permute(k(:,1:4:end,:),[2 1 3]);
Fin{2} = zeros([60 96 32]); Fin{2}(4:4:end,:,:) = permute(k(:,3:4:end,:),[2 1 3]);
Fdpg2 = dpg_recon( Fin, N2, 2, 2);

