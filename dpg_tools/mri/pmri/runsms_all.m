if ~exist('verbose'), verbose=1; end;


if ~exist('fname');
  determine_fname;
end;

if ~exist('hdr'); hdr = parse_measdat_hdr2( fname ); end;
R = hdr.R;
slcindx = hdr.slcindx;

if ~exist('prot')
  prot = eval_ascconv( fname );
end;
evp.NSlcMeas = length(slcindx);
% [SMS, FOVshift, NSlc, Ngroup, PhaseShift, SliceSep] = mrir_array_SMS_recon_params(prot,evp);
[SMS, FOVshift, NSlc, NcSlc, PhaseShift] = sms_recon_parms(hdr);
R = hdr.R;
% NcSlc = NSlc / SMS;                     % number of collapsed-slice measurements

[~,idx] = sort(slcindx);
slices = 1:length(slcindx);
indx = reshape( slices(idx), [ length(slices)/SMS SMS ])
clear idx slices

[~,idx] = sort(min(indx'));
indx = indx(idx,:)

if ~exist('dpgkernel')
  dpgkernel = '3x5';
end;


if ~exist('ksrc','var')
  % load the measured data mapping
  if iscell(k)
    ksrc = k{length(k)};
  else
    ksrc = k;
  end;
  if isempty(v)
    gg = mrir_regrid_trapezoid_scottdata( ksrc.image(:,coils,:,1,:,1,1,1,Tcnt,:), meas.prot, regrid_trapezoid_prep);
  else
    gg = tmult(ksrc.image(:,coils,:,1,:,1,1,1,Tcnt,:),v',1);
  end;
  if (size(gg, length( size(gg) ))==2)
    gg = sum( gg, length(size(gg)) );
  end;
end;

if ( exist('runmod') && ( runmod == 1 ) )
    % data on 2016_06_14 had a corrupted channel-13.  So, remove it from
    % processing here.
    coils = [1:12 14:double(ksrc.image.dataSize(2))];
else
  runmod = 0;
  coils = [1:double(ksrc.image.dataSize(2))];
end;

% collapsed data deinterleaving
%%% this is a bit confusing (and a source or a bug in May 2016) ...
%%% the interleaved slice ordering is determined by the total number of slices.  
%%% but the number of sampled slices is actually the collapsed NcSlc number.
if ( floor(NSlc/2) == (NSlc/2) )
  slcorder = [ 2:2:NcSlc 1:2:NcSlc ];
else
  slcorder = [ 1:2:NcSlc 2:2:NcSlc ];
end;


for Tcnt = 1:ksrc.image.dataSize(9);

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

  Fs1 = zeros(sz([2 1 3 4]));
  for slc=1:NSlc/SMS;

    Fin.p = zeros(sz([2 1 4]));
    Fin.p(ggi(1:2:end),:,:) = permute( squeeze(k_data_deblur(:,1:2:end,slc,:)),[2 1 3]);
    Fin.n = zeros(sz([2 1 4]));
    Fin.n(ggi(2:2:end),:,:) = permute( squeeze(k_data_deblur(:,2:2:end,slc,:)),[2 1 3]);
    
    for cnt=1:SMS;
      rslc = (cnt-1)*(NSlc/SMS) + slc;
      Fs1(:,:,rslc,:) = dpg_recon( Fin, Nsms_lb{slc}{cnt}, R, dpgkernel );
      % remove the relative phase shift
      Fs1(:,:,rslc,:) = tmult( Fs1(:,:,rslc,:), diag(exp(j*-(cnt-1)*PhaseShift* ...
                                                        (0:(size(Fs1,1)-1))/R )),1);

      if (verbose)
        im( sqrt(sum(abs( fif( Fs1(:,:,rslc,:) ) ).^2,4)), gcf );
        pause(0.5);
      end;
    end;
    
  end;

  Fs1 = permute(Fs1,[2 1 3 4]);
  m = max(sz([1 2]));
  tmp = zeros([ m m sz([3 4]) ]);
  tmp( :, floor((m-size(Fs1,2))/2)+(1:size(Fs1,2)),:,:) = Fs1;
  Is2 = sqrt(sum(abs( flipdim(fif(tmp),1) ).^2,4));
  
  % use 'Is3' to denote that this is using the C3 ACS data
  savend( sprintf('Is2_%03d.nd',Tcnt), Is2([end 1:(end-1)],:,:), 'flt');

  if (verbose)
    im( cascade( Is2 ), gcf );
    keyboard;
  end;
end;

