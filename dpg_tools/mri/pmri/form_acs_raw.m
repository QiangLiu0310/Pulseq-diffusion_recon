%% generate an ACS image from the first/standard ACS data set.

if ~exist('hdr','var');
  hdr = parse_measdat_hdr(fname); 
  hdr.NRefLin = k.refscan.dataSize(3);
  hdr.NColMeas = k.refscan.dataSize(1);
end;

v = extract_vrgf(fname);

fid = fopen(fname,'r');

segln = hdr.NRefLin / hdr.R;

minln = Inf;
for cnt=1:hdr.NSlc*hdr.R;   %only work for R=3 or R>3, because R=2 the ACS is single shot

  mpPC  = k.refscanPC.memPos( (cnt-1)*3 + 1 );

  if strcmp( hdr.data_version, 'VB' )
    sz2 = [ ( hdr.NColMeas + 16) ...
            hdr.ncoils ...
            3 + segln ];
    nEl = prod(sz2);
  else
    sz2 = [ hdr.NColMeas ... % sMDH.ushSamplesInScan ...
            hdr.ncoils ...
            3 + hdr.NRefLin/hdr.R ... % ascconv.sKSpace.lPhaseEncodingLines ...
            ... % hdr.NSlc ... % ascconv.sSliceArray.lSize ...
            ];
    mdhsz = 192 / 8;                      % number of 'float' elements in the MDH to discard
    cdhsz = 32 / 8;                       % number of 'float' elements in the ChanHdr to discard
    
    nEl =  (mdhsz + ( cdhsz+sz2(1) )*sz2(2))*prod(sz2(3:end));
  end;

  fseek(fid,mpPC,'bof');
  a = fread(fid,nEl*2,'float=>real*4');

  if strcmp( hdr.data_version, 'VB' )
    b = reshape( complex(a(1:2:nEl*2),a(2:2:nEl*2)), sz2 );
    b = b(17:end,:,:,:,:);
    
  else
    b = reshape( complex(a(1:2:end),a(2:2:end)), [ (mdhsz + ( cdhsz+sz2(1) )*sz2(2)) sz2(3:end) ]);
  
    bsz = size(b); bsz(1) = bsz(1) - mdhsz;
    b2 = reshape( b( (mdhsz+1):end,:), bsz ); % trim off the MDH values

    b3 = reshape( b2, sz2 + [ cdhsz zeros([1 length(sz2)-1]) ] );

    b = reshape( b3( (cdhsz+1):end,:), sz2 );
  end;

  b(:,:,1:2:end,:) = flipdim( b(:,:,1:2:end,:),1);
  if (v ~= 1), b = tmult( b, v.', 1); end;
    
  kdat = permute(b,[ 1 3 4 2 5]);
  vdat = kdat(:,1:3,:,:,:);
  kdat = kdat(:,4:end,:,:,:);

  mpRef = k.refscan.memPos( (cnt-1)*segln + 1 );  

  fseek(fid,mpRef,'bof');
  if strcmp( hdr.data_version, 'VB' )
    sMDH = ice_read_mdh_vb17(fid);
  else
    sMDH = ice_read_mdh_vd13(fid);
  end;
  fprintf(' Ref: %d %d %d \n',sMDH.sLC.ushSlice, sMDH.sLC.ushLine, sMDH.sLC.ushSeg );

  if ( (sMDH.sLC.ushLine + 1) < minln)
    minln = sMDH.sLC.ushLine + 1;
  end;

  Sk0 = squeeze(mean( vdat(:,[1 3],:,:),2 ));
  Sk1 = squeeze(mean( vdat(:,[ 2 ],:,:),2 ));

  Sk0 = ffts(ifft(ffts(Sk0,1),[],1),1);
  Sk1 = ffts(ifft(ffts(Sk1,1),[],1),1);

  phz = comp_local_pc( Sk1, Sk0 );
  kdat2 = phzapply( permute(kdat,[2 1 4 3]), phz );

  Craw( 1 + sMDH.sLC.ushLine + (0:(segln-1))*hdr.R, :, sMDH.sLC.ushSlice + 1, : ) = kdat2;

end;

Craw = Craw( minln:end,:,:,:) ;
fclose(fid);

if isfield(hdr,'SMS')
  [SMS, FOVshift, NSlc, Ngroup, PhaseShift] = sms_recon_parms(hdr);
 
  [~,idx] = sort(hdr.slcindx);
  slc = 1:length(hdr.slcindx);
  indx = reshape( slc(idx), [ length(slc)/SMS SMS ])
  clear idx slc
 
  if (SMS>1)
    [~,idx] = sort(min(indx'));
  else
    [~,idx] = sort(indx);
  end;
  indx = indx(idx,:)
end;
