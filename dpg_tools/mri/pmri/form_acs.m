%determine_fname
if ~exist('k','var')
  k = mapVBVD(fname);
end;

if (1)
  [meas.prot, meas.evp] = read_meas_prot(fname);
  regrid_trapezoid_prep = mrir_regrid_trapezoid_prep(meas.prot, size(k.refscan(:,:,:,1, 01,1,1,1,1,1,1), 1));
  v = [];
else
  v = extract_vrgf(fname);
end;

hdr = parse_measdat_hdr(fname);

if ~exist('prot')
  prot = eval_ascconv( fname );
end;
evp.NSlcMeas = length(hdr.slcindx);
% [SMS, FOVshift, NSlc, Ngroup, PhaseShift, SliceSep] = mrir_array_SMS_recon_params(prot,evp);
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

if iscell(k)
  ksrc = k{2};
else
  ksrc = k;
end;
if isempty(v)
  A = mrir_regrid_trapezoid_scottdata( ksrc.refscan(:,:,:,1,:,1,1,1,1,1,1), meas.prot, regrid_trapezoid_prep);
  B = mrir_regrid_trapezoid_scottdata( ksrc.refscan(:,:,:,1,:,1,1,1,1,1,2), meas.prot, regrid_trapezoid_prep);
else
  A = tmult( ksrc.refscan(:,:,:,1,:,1,1,1,1,1,1),v',1);
  B = tmult( ksrc.refscan(:,:,:,1,:,1,1,1,1,1,2),v',1);
end;
A = permute(A,[1 3 5 2 4]);
B = permute(B,[1 3 5 2 4]);

C = ifi( permute( phzshift( permute(fif(A),[2 1 3 4]), permute(fif(B),[2 1 3 4]),'nofft' ), [2 1 3 4]) );

C2(:,:,hdr.slcindx,:) = C;
A2(:,:,hdr.slcindx,:) = A;
B2(:,:,hdr.slcindx,:) = B;
