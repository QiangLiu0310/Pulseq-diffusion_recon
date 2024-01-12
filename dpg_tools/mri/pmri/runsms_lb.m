
%%% run DPG SMS with leak-block...

%% Inputs
hdr = parse_measdat_hdr2( fname );

if ~exist('prot')
  prot = eval_ascconv( fname );
end;
evp.NSlcMeas = length(hdr.slcindx);

% [SMS, FOVshift, NSlc, Ngroup, PhaseShift, SliceSep] = mrir_array_SMS_recon_params(prot,evp);
[SMS, FOVshift, NSlc, NcSlc, PhaseShift] = sms_recon_parms(hdr);
R = hdr.R;
% NcSlc = NSlc / SMS;                     % number of collapsed-slice measurements


if ~exist('verbose','var'), verbose = 0; end;
if ~exist('dbg','var'), dbg = 0; end;
if ~exist('solver','var'), solver = 'cgsr'; end;
if ~exist('n','var'), n = 0; end;
if ~exist('chi','var'), chi = 1e-6; end;
if ~exist('v','var'), v = []; end;

if exist('C3','var')
  sz = size(C3);
  k_acs1 = zeros(sz); k_acs1(:,(1:sz(2)),:,:) = A3;
  k_acs2 = zeros(sz); k_acs2(:,(1:sz(2)),:,:) = B3;
  k_trgt = zeros(sz); k_trgt(:,(1:sz(2)),:,:) = C3;
else
  sz = size(C2);
  k_acs1 = A2;
  k_acs2 = B2;
  k_trgt = C2;
end;

phzabs = exp(j*(PhaseShift*(0:(size(k_trgt,2)-1))/R - (0.0*pi) - PhaseShift ) );

% FOV shifting and phase modulating...
sz = size(k_acs1);

for cnt=0:(SMS-1),
  k_acs1(:,:,cnt*NcSlc+(1:NcSlc),:) = tmult( k_acs1(:,:,(cnt*NcSlc)+(1:NcSlc),:), diag(phzabs.^cnt), 2);

%   k_acs1_collapsed = k_acs1_collapsed + k_acs1(:,:,cnt*NcSlc+(1:NcSlc),:);
end;


% k_acs2_collapsed = zeros([ sz(1:2) NcSlc sz(4)]);
for cnt=0:(SMS-1),
    k_acs2(:,:,cnt*NcSlc+(1:NcSlc),:) = tmult( k_acs2(:,:,(cnt*NcSlc)+(1:NcSlc),:), diag(phzabs.^cnt), 2);
    
%     k_acs2_collapsed = k_acs2_collapsed + k_acs2(:,:,cnt*NcSlc+(1:NcSlc),:);
end;

if (dbg)
  k_acs1_collapsed = zeros([ sz(1:2) NcSlc sz(4)]);
  k_acs2_collapsed = zeros([ sz(1:2) NcSlc sz(4)]);
  for cnt=0:(SMS-1),
    k_acs1_collapsed = k_acs1_collapsed + k_acs1(:,:,cnt*NcSlc+(1:NcSlc),:);
    k_acs2_collapsed = k_acs2_collapsed + k_acs2(:,:,cnt*NcSlc+(1:NcSlc),:);
  end;
end;

% k_acs2(:,:,1:3,:) = CaipiShift_K( k_acs2(:,:,1:3,:),  0, PhaseShift );
% k_acs2(:,:,4:6,:) = CaipiShift_K( k_acs2(:,:,4:6,:), -1, PhaseShift );
% k_acs2(:,:,7:9,:) = CaipiShift_K( k_acs2(:,:,7:9,:), -2, PhaseShift );
% 
% k_acs2_collapsed = k_acs2(:,:,1:3,:) + ...
%     tmult( k_acs2(:,:,4:6,:), diag(phzabs), 2) + ...
%     tmult( k_acs2(:,:,7:9,:), diag(phzabs.^2), 2);

for cnt=0:(SMS-1),
  k_trgt(:,:,(cnt*NcSlc)+(1:NcSlc),:) = tmult( k_trgt(:,:,(cnt*NcSlc)+(1:NcSlc),:), diag(phzabs.^cnt), 2);
end;

% collapsed data deinterleaving
if ( floor(NSlc/2) == (NSlc/2) )
  slcorder = [ 2:2:NcSlc 1:2:NcSlc ];
else
  slcorder = [ 1:2:NcSlc 2:2:NcSlc ];
end;

if exist('kraw','var')
  k_data = zeros([size(kraw,1) size(kraw,2)-3 NcSlc size(kraw,4)]);
  k_data(:,1:end,slcorder,:) = kraw(:,4:end,2*NSlc+(1:NcSlc),:); 


  % collapsed data Deblurring
  k_data_deblur = zeros(size(k_data));
  k_data_deblur(:,:,:,:) = CaipirinhaDeblur_v3_wsh( k_data(:,:,:,:), prot, evp );
% for SlcCount = 1:size(k_data,3)
%   k_data_deblur(:,:,SlcCount,:) = CaipiShift_K(k_data(:,:,SlcCount,:), ...
%                                                (SlcCount-1-IsoSlice)/5,-PhaseShift);
%   % the "-1" is to match counting from 'zero' in the sequence
% end
elseif exist('k','var')
  Tcnt = 2;
  if iscell(k)
    ksrc = k{length(k)};
  else
    ksrc = k;
  end;
  if ~exist('coils','var'); coils = 1:ksrc.image.dataSize(2); end;
  if isempty(v)
    gg = mrir_regrid_trapezoid_scottdata( ksrc.image(:,coils,:,1,:,1,1,1,Tcnt,:), meas.prot, regrid_trapezoid_prep);
  else
    gg = tmult(ksrc.image(:,coils,:,1,:,1,1,1,Tcnt,:),v',1);
  end;
  if (size(gg, length( size(gg) ))==2)
    gg = sum( gg, length(size(gg)) );
  end;

  ggi = find( sum(sum( gg(:,:,:,1),1) ,2) );

  sz = size(gg);
  sz = sz([1 3 5 2]);                   % [kx ky z coils]

  k_data = permute(gg(:,:,ggi,1,indx(:,1)),[1 3 5 2 4]);
  k_data(:,:,slcorder,:) = k_data;
  if isfield(prot,'sSliceAcceleration')
    % for data collected from the 'smsprod' seq, (eg 7T Terra data), data is already 'deblurred'
    k_data_deblur = k_data;
  else
    % this line is needed for older 'mgh' versions of the sequence
    k_data_deblur = CaipirinhaDeblur_v3_wsh( k_data(:,:,:,:), prot, evp );
  end;
else
  k_data_deblur = [];
end;

if (dbg)
  % tmp = zeros([ sz(1) floor((sz(2)+1)/R) NcSlc sz(4)]);
  % tmp(:,1:2:end,:,:) = k_acs1_collapsed(:,2:2*R:(size(tmp,2)*R),:,:);
  % tmp(:,2:2:end,:,:) = k_acs2_collapsed(:,(2+R):2*R:(size(tmp,2)*R),:,:);
  tmp = zeros(size(k_acs1_collapsed));
  tmp(:,1:2:end,:,:) = k_acs1_collapsed(:,1:2:end,:,:);
  tmp(:,2:2:end,:,:) = k_acs1_collapsed(:,2:2:end,:,:);
  
  % im([ (tmp) (k_data_deblur) ],gcf);  pause(0.1);
  cmprindx = floor((size(k_data_deblur,2)-size(tmp,2)/hdr.R)/2)+(1:size(tmp,2)/hdr.R);
  im([ (tmp(:,1:hdr.R:end,:,:)) (k_data_deblur(:,cmprindx,:,:,:)) ],gcf)
  
  % im([ fif(k_acs1_collapsed) fif(k_data_deblur) ],gcf); pause(0.5);
  keyboard;
end;

if ~exist('dpgkernel')
  dpgkernel = '3x5';
end;

disp(sprintf([' using solver: "%s", normalize: %d, chi: %e'],solver,n,chi));

clear Nsms_lb;
fprintf('%d',mod(1:NcSlc,10)); fprintf('\n');
for slc = 1:NcSlc,
  fprintf('.');

  if ( ~exist('Nsms_lb') || (length(Nsms_lb)<slc) || isempty(Nsms_lb{slc}) )
    kin.p = permute( k_acs1(:,:,slc:NcSlc:end,:),[ 2 1 3 4]);
    kin.n = permute( k_acs2(:,:,slc:NcSlc:end,:),[ 2 1 3 4]);
    kin.target = permute( k_trgt(:,:,slc:NcSlc:end,:),[ 2 1 3 4]);
    [~,~,Nsms_lb{slc}] = dpg_sms( size(kin.p), kin, vec(1:size(kin.p,1)), 'kernel',dpgkernel,'dks',R,'solver',solver,'normalize',n,'chi',chi,'fast'); % ,'report',1);
  end;

  if ~isempty(k_data_deblur)
    Fin.p = zeros([ R*2*size(k_data_deblur,2)/2 sz([1 4]) ]);
    Fin.p(1:2*R:end,:,:) = permute( squeeze(k_data_deblur(:,1:2:end,slc,:)),[2 1 3]);
    Fin.n = zeros([ R*2*size(k_data_deblur,2)/2 sz([1 4]) ]);
    Fin.n((1+R):2*R:end,:,:) = permute( squeeze(k_data_deblur(:,2:2:end,slc,:)),[2 1 3]);

    for cnt=1:SMS;
      rslc = (cnt-1)*NcSlc + slc;
      Fsms2(:,:,rslc,:) = dpg_recon( Fin, Nsms_lb{slc}{cnt}, R, dpgkernel );
    end;
    if (verbose)
      im( cascade(sqrt(sum(abs( fif( Fsms2(:,:,(0:SMS-1)*NcSlc+slc,:) ) ).^2,4)),[1 SMS]), gcf );
      pause(0.5);
    end;
  end;

end;
fprintf('\n');

% shift the SMS slices back:
% Fsms(:,:,4:6,:) = permute( CaipiShift_K( permute(Fsms(:,:,4:6,:),[2 1 3 4]), 1, PhaseShift ),[2 1 3 4]);
% Fsms(:,:,7:9,:) = permute( CaipiShift_K( permute(Fsms(:,:,7:9,:),[2 1 3 4]), 2, PhaseShift ),[2 1 3 4]);

phzabs = exp(j*(PhaseShift*(0:(size(Fsms2,1)-1))/R - (0.0*pi) - PhaseShift ) );
for cnt=0:(SMS-1),
  Fsms2(:,:,(cnt*NcSlc)+(1:NcSlc),:) =  tmult( Fsms2(:,:,(cnt*NcSlc)+(1:NcSlc),:), diag(conj(phzabs).^cnt), 1);
end;

Isms2 = sqrt(sum(abs( fif(Fsms2) ).^2,4));
