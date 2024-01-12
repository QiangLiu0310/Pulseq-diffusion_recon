%% should be run after run_2kernal_lb.m

if ~exist('verbose'), verbose=0; end;

if ~exist('fname');
  determine_fname;
end;

if ~exist('hdr')
  hdr = parse_measdat_hdr( fname );
end;
R = hdr.R;
slcindx = hdr.slcindx;

if ~exist('prot')
  prot = eval_ascconv( fname );
end;
evp.NSlcMeas = length(slcindx);
% [SMS, FOVshift, NSlc, Ngroup, PhaseShift, SliceSep] = mrir_array_SMS_recon_params(prot,evp);
[SMS, FOVshift, NSlc, NcSlc, PhaseShift] = sms_recon_parms(hdr);
NcSlc = NSlc / SMS;                     % number of collapsed-slice measurements

[~,idx] = sort(slcindx);
slices = 1:length(slcindx);
indx = reshape( slices(idx), [ length(slices)/SMS SMS ])
clear idx slices

[~,idx] = sort(min(indx'));
indx = indx(idx,:)

if ~exist('k')
  % load the measured data mapping
  k = mapVBVD(fname);
  % k.image.flagIgnoreSeg = 1;

  v = extract_vrgf( fname );
end;

% collapsed data deinterleaving
if ( floor(NSlc/2) == (NSlc/2) )
  slcorder = [ 2:2:NcSlc 1:2:NcSlc ];
else
  slcorder = [ 1:2:NcSlc 2:2:NcSlc ];
end;

if ~exist('Np')
  Np = load_gparms('Np.nd');
end;

z = [];
for cnt=1:length(Np{1})
  dk = diff(find( Np{1}(cnt).pattern == '*' ));
  z(dk) =  1;
end;

if ( z( 2*hdr.R ) == 0 )
  disp(['for this script to run properly, GRAPPA parameters in Np need to cover both R and 2R accelerations']);
  error(['need to regenerate Np parameters'])
end;

if ~exist('w','var')
  load w_2kernel_coef
end;

if iscell(k)
  ksrc = k{length(k)};
else
  ksrc = k;
end;
if ( exist('runmod') && ( runmod == 1 ) )
    % data on 2016_06_14 had a corrupted channel-13.  So, remove it from
    % processing here.
    coils = [1:12 14:double(ksrc.image.dataSize(2))];
else
  runmod = 0;
  coils = [1:double(ksrc.image.dataSize(2))];
end;

nky = size(k_data_deblur,2);

for Tcnt = 1:ksrc.image.dataSize(9);
  fprintf(':');

  if isempty(v)
    gg = mrir_regrid_trapezoid_scottdata( ksrc.image(:,coils,:,1,:,1,1,1,Tcnt,:), meas.prot, regrid_trapezoid_prep);
    vv = mrir_regrid_trapezoid_scottdata( ksrc.phasecor(:,coils,:,1,:,1,1,1,Tcnt,:), meas.prot, regrid_trapezoid_prep);
  else
    gg = tmult(ksrc.image(:,coils,:,1,:,:,1,1,Tcnt,:),v',1);
    % vv = tmult(ksrc.phasecor(:,coils,:,1,:,:,1,1,Tcnt,:),v',1);
    vv = tmult( ksrc.phasecor(:,coils,:,1,indx(:,1),:,1,1,Tcnt,:), v', 1 );
  end;
  if (size(gg, length( size(gg) ))==2)
    gg = sum( gg, length(size(gg)) );
  end;
  ggi = find( sum(sum( gg(:,:,:,1),1) ,2) );

  sz = size(gg);
  sz = sz([1 3 5 2]);                   % [kx ky z coils]

  k_data = permute(gg(:,:,ggi,1,indx(:,1)),[1 3 5 2 4]);
  k_data(:,:,slcorder,:) = k_data;

  % v_data = kraw(:,1:3,2*NSlc+(1:NcSlc),:); 
  v_data = permute( squeeze( vv(:,:,end,1,:,:,:)),[1 4 3 2 5]);
  S1 = sum(v_data(:,:,:,:,1),2)./sum(v_data(:,:,:,:,1)~=0,2);
  S2 = sum(v_data(:,:,:,:,2),2)./sum(v_data(:,:,:,:,2)~=0,2);
  S1 = ffts(ifft(ffts( S1,1),[],1),1);
  S2 = ffts(ifft(ffts( S2,1),[],1),1);

  if isfield(prot,'sSliceAcceleration')
    % for data collected from the 'smsprod' seq, (eg 7T Terra data), data is already 'deblurred'
    k_data_deblur = k_data;
  else
    % this line is needed for older 'mgh' versions of the sequence
    k_data_deblur = CaipirinhaDeblur_v3_wsh( k_data(:,:,:,:), prot, evp );
  end;

  phz = comp_local_pc( S1, S2 );
  k_data_gc = permute( phzapply( permute( k_data,[2 1 3 4]), phz ), [2 1 3 4] );

  if isfield(prot,'sSliceAcceleration')
    % for data collected from the 'smsprod' seq, (eg 7T Terra data), data is already 'deblurred'
    k_data_gc_deblur = k_data_gc;
  else
    % collapsed data Deblurring
    k_data_gc_deblur = CaipirinhaDeblur_v3_wsh( k_data_gc(:,:,:,:), prot, evp );
  end;
  
  for slc = 1:NcSlc;
    fprintf('-');
    in1 = zeros([ sz(1) nky*hdr.R sz(4)]);

    in1(:,hdr.R*(0:size(k_data_gc_deblur,2)-1)+1,:) = squeeze( k_data_gc_deblur(:,:,slc,:) );
    
    iii = 1:size(k_trgt,2);
    for cnt=1:SMS,
      in2{cnt} = zeros( sz([1 1 4]) );
      in2{cnt}(:,(sz(1)-length(iii))/2+iii,:) = squeeze( k_trgt(:,:,slc+(cnt-1)*NcSlc,:) ); 
    end;
    % [~,w{slc}] =  MultisliceGRAPPA_2kernal_leakBlock( in1, in2, [ 7 5 1 R ], 'full', prot);

    tmp = squeeze(MultisliceGRAPPA_2kernal_leakBlock( in1, w{slc}, [ 5 3 1 hdr.R ]) );

    % phzabs = exp(-j*(PhaseShift*(0:(size(tmp,2)-1))/R - (0.0*pi) - PhaseShift) );
    % for cnt=1:SMS,
    %   tmp(:,:,:,cnt) = tmult( tmp(:,:,:,cnt), diag( phzabs.^(cnt-1)), 2);
    % end;

    phzabs = exp(j*(PhaseShift*(0:(size(tmp,2)-1))/hdr.R - (0.0*pi) - PhaseShift ) );
    slcgrp = slc + [ 0:(SMS-1) ]*NcSlc;

    sz_in1 = size(in1);
    for cnt=1:SMS,
      curslc = slcgrp(cnt);
      Fa1 = recongrappa_multik(sz_in1([2 1 3]),permute(tmp(:,1:2*hdr.R:end,:,cnt),[2 1 3]),1:2*hdr.R:sz_in1(2),'kernel','2x5','N',Np{curslc});
      Fb1 = recongrappa_multik(sz_in1([2 1 3]),permute(tmp(:,(1+hdr.R):2*hdr.R:end,:,cnt),[2 1 3]),(1+hdr.R):2*hdr.R:sz_in1(2),'kernel','2x5','N',Np{curslc});
    
      Fa2 = permute(tunfold(fif(Fa1),2),[2 1]);
      Fb2 = permute(tunfold(fif(Fb1),2),[2 1]);

      Fb3 = phzshift( Fa2, Fb2,{'nofft','nocombo'} );
      Fb4 = ifi( trefold(permute(Fb3,[2 1]),sz_in1([2 1 3]),2) );

      Fc1 = zeros(size(Fa1));
      Fc1( 1:2*hdr.R:end, :, : ) = Fa1( 1:2*hdr.R:end, :, : );
      Fc1( (1+hdr.R):2*hdr.R:end, :, : ) = Fb4( (1+hdr.R):2*hdr.R:end, :, : );

      Fd1 = recongrappa_multik(size(Fc1),Fc1,[],'kernel','2x5','dks',hdr.R,'N',Np{curslc});

      F2klb(:,:,curslc,:) = Fd1;
      F2klb(:,:,curslc,:) =  tmult( F2klb(:,:,curslc,:), diag(conj(phzabs).^(cnt-1)), 1);
    end;
  end;

  I2k = zeros([ size(F2klb,2)*[1 1] size(F2klb,3) size(F2klb,4) ]);
  I2k( (size(I2k,1)-size(F2klb,1))/2+(1:size(F2klb,1)),:,:,:) = flipdim(F2klb,2);
  I2k = sqrt(sum(abs( fif(I2k) ).^2,4));

  fprintf('\n');
  savend(sprintf('I2k_%03d.nd',Tcnt),I2k,'flt');
  if (verbose); keyboard; end;
end;
