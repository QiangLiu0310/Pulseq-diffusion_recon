%% should be run after runsms_lb.m
if ~exist('verbose','var'); verbose = 0; end;
if ~exist('NcSlc','var'); NcSlc = hdr.NSlc / hdr.SMS; end;
if exist('C3')
  k_trgt = C3;
else
  k_trgt = C2;
end;
sz = size(k_trgt);

if ~exist('slices','var'); slices = 1:hdr.NSlc; end;
if ~exist('chi','var'); chi = 1e-6; end;
if ~exist('eta','var'); eta = 1; end;

if ~exist('prot')
  prot = eval_ascconv( fname );
end;
evp.NSlcMeas = length(hdr.slcindx);

[SMS, FOVshift, NSlc, NcSlc, PhaseShift] = sms_recon_parms(hdr);
R = hdr.R;

phzabs = exp(j*(PhaseShift*(0:(size(k_trgt,2)-1))/hdr.R - (0.0*pi) - PhaseShift ) );

for cnt=0:(SMS-1),
  k_trgt(:,:,(cnt*NcSlc)+(1:NcSlc),:) = tmult( k_trgt(:,:,(cnt*NcSlc)+(1:NcSlc),:), diag(phzabs.^cnt), 2);
end;


% collapsed data deinterleaving
if ( floor(hdr.NSlc/2) == (hdr.NSlc/2) )
  slcorder = [ 2:2:NcSlc 1:2:NcSlc ];
else
  slcorder = [ 1:2:NcSlc 2:2:NcSlc ];
end;

if ~exist('Np','var')
  if exist('Np.nd','file')
    Np = load_gparms('Np.nd');
  else
    fprintf('%d',mod(1:length(slices),10));fprintf('\n');
    for slc = slices; % 1:NSlc;
      fprintf('o');
      [~,~,Np{slc}]=recongrappa_multik(sz([2 1 4]),permute(k_trgt(:,:,slc,:),[2 1 4 3]),[],'kernel','2x5','dks',hdr.R*[1 2],...
                                       'chi',chi,'eta',eta );
    end;
    fprintf('\n');
  end;
end;

z = [];
for cnt=1:length(Np{slices(1)})
  dk = diff(find( Np{slices(1)}(cnt).pattern == '*' ));
  z(dk) =  1;
end;

if ( (length(z)< 2*hdr.R) || ( z( 2*hdr.R ) == 0 ) )
  disp(['for this script to run properly, GRAPPA parameters in Np need to cover both R and 2R accelerations']);
  error(['need to regenerate Np parameters'])
end;

if exist('kraw','var')
  k_data = zeros([size(kraw,1) size(kraw,2)-3 NcSlc size(kraw,4)]);
  % k_data(:,1:end,slcorder,:) = kraw(:,4:end,2*NSlc+(1:NcSlc),:); 
  k_data(:,1:end,slcorder,:) = kraw(:,4:end,(1:NcSlc),:); 

  % v_data = kraw(:,1:3,2*NSlc+(1:NcSlc),:); 
  v_data = zeros([size(kraw,1) 3 NcSlc size(kraw,4)]);
  v_data(:,:,slcorder,:) = kraw(:,1:3,(1:NcSlc),:); 

  S1 = mean(v_data(:,[1 3],:,:),2); 
  S2 = mean(v_data(:,[ 2 ],:,:),2);
  S1 = ffts(ifft(ffts( S1,1),[],1),1);
  S2 = ffts(ifft(ffts( S2,1),[],1),1);

else
  if (k.image.dataSize(9) > 1),
    Tcnt = 1;
  else
    Tcnt = 1;
  end;

  if iscell(k)
    ksrc = k{length(k)};
  else
    ksrc = k;
  end;

  if ~exist('coils','var'); coils = 1:ksrc.image.dataSize(2); end;
  if isempty(v)
    gg = mrir_regrid_trapezoid_scottdata( ksrc.image(:,coils,:,1,:,1,1,1,Tcnt,:), meas.prot, regrid_trapezoid_prep);
    vv = mrir_regrid_trapezoid_scottdata( ksrc.phasecor(:,coils,:,1,:,:,1,1,Tcnt,:), meas.prot, regrid_trapezoid_prep);
  else
    gg = tmult(ksrc.image(:,coils,:,1,:,:,1,1,Tcnt,:),v',1);
    % vv = tmult(ksrc.phasecor(:,coils,:,1,:,:,1,1,Tcnt,:),v',1);
    vv = tmult( ksrc.phasecor(:,coils,:,1,:,:,1,1,Tcnt,:), v', 1 );
  end;
  if (size(gg, length( size(gg) ))==2)
    gg = sum( gg, length(size(gg)) );
  end;
  ggi = find( sum(sum( gg(:,:,:,1),1) ,2) );

  sz = size(gg);
  sz = sz([1 3 5 2]);                   % [kx ky z coils]

  v_data = permute( squeeze( vv(:,:,end,1,indx(:,1),:,:)),[1 4 3 2 5]);
  S1 = sum(v_data(:,:,:,:,1),2)./sum(v_data(:,:,:,:,1)~=0,2);
  S2 = sum(v_data(:,:,:,:,2),2)./sum(v_data(:,:,:,:,2)~=0,2);
  S1 = ffts(ifft(ffts( S1,1),[],1),1);
  S2 = ffts(ifft(ffts( S2,1),[],1),1);


  k_data = permute(gg(:,:,ggi,1,indx(:,1)),[1 3 5 2 4]);
  k_data(:,:,slcorder,:) = k_data;
  if isfield(prot,'sSliceAcceleration')
    % for data collected from the 'smsprod' seq, (eg 7T Terra data), data is already 'deblurred'
    k_data_deblur = k_data;
  else
    % this line is needed for older 'mgh' versions of the sequence
    k_data_deblur = CaipirinhaDeblur_v3_wsh( k_data(:,:,:,:), prot, evp );
  end;
  
end;


phz = comp_local_pc( S1, S2 );

k_data_gc = permute( phzapply( permute( k_data,[2 1 3 4]), phz ), [2 1 3 4] );

if isfield(prot,'sSliceAcceleration')
  % for data collected from the 'smsprod' seq, (eg 7T Terra data), data is already 'deblurred'
  k_data_gc_deblur = k_data_gc;
else
  % collapsed data Deblurring
  % this is needed for older 'mgh' versions of the sequence
  k_data_gc_deblur = CaipirinhaDeblur_v3_wsh( k_data_gc(:,:,:,:), prot, evp );
end;

if exist('runmod')
  if ( runmod == 1 )
    k_data_gc_deblur = k_data_gc_deblur(:,:,:,coils);
  end;
end;

nky = size(k_data_deblur,2);

fprintf('%d',mod(1:length(slices),10)); fprintf('\n');
for slc = slices(1:end/hdr.SMS); % 1:NcSlc;
  fprintf('.');

  in1 = zeros([ size(k_data_gc_deblur,1) nky*hdr.R sz(4)]);

  in1(:,hdr.R*(0:size(k_data_gc_deblur,2)-1)+1,:) = squeeze( k_data_gc_deblur(:,:,slc,:) );

  if (~exist('w','var') || (length(w)<slc) || isempty(w{slc}) ),
    iii = 1:size(k_trgt,2);
    sz = size(k_trgt);
    for cnt=1:hdr.SMS,
      in2{cnt} = zeros( sz([1 2 4]) );
      in2{cnt}(:,(sz(2)-length(iii))/2+iii,:) = squeeze( k_trgt(:,:,slc+(cnt-1)*NcSlc,:) ); 
    end;
    [~,w{slc}] =  MultisliceGRAPPA_2kernal_leakBlock( in2{1}, in2, [ 5 3 1 hdr.R ], 'full', prot);
  end;

  tmp = squeeze(MultisliceGRAPPA_2kernal_leakBlock( in1, w{slc}, [ 5 3 1 hdr.R ]) );


 phzabs = exp(j*(PhaseShift*(0:(size(tmp,2)-1))/hdr.R - (0.0*pi) - PhaseShift) );
 for cnt=1:hdr.SMS,
   tmp(:,:,:,cnt) = tmult( tmp(:,:,:,cnt), diag( phzabs.^(cnt-1)), 2);
 end;

  phzabs = exp(j*(PhaseShift*(0:(size(tmp,2)-1))/hdr.R - (0.0*pi) - PhaseShift ) );
  slcgrp = slc + [ 0:(hdr.SMS-1) ]*NcSlc;

  sz_in1 = size(in1);
  for cnt=1:hdr.SMS,
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

  if verbose, keyboard; end;
end;
fprintf('\n');

% phzabs = exp(j*(PhaseShift*(0:(size(F2klb,1)-1))/hdr.R - (0.0*pi) - PhaseShift ) );
% for cnt=0:(SMS-1),
%   F2klb(:,:,(cnt*NcSlc)+(1:NcSlc),:) =  tmult( F2klb(:,:,(cnt*NcSlc)+(1:NcSlc),:), diag(conj(phzabs).^cnt), 1);
% end;
