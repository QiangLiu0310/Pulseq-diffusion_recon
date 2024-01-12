function a3 = rfmt_ima( fname, slc, m_in, n_in );
% converts a mosaic-ed DICOM image into a 3D image volume
%
% fname - the input DICOM file (or mosaic'd image)
% slc   - the number of (true) slices in the image

if ischar( fname )
  if isdicom( fname )
    a = dicomread( fname );
  else
    fid = fopen(fname);
    a = fread(fid,inf,'short');
    fclose(fid);
    
    m = sqrt(length(a));
    
    a = reshape(a,[m m]).';
  end;
else
  a = fname;
end;


[m n] = size(a);

if (nargin==1)
  if isdicom( fname ),
    hdr = dicominfo(fname);
    if isfield(hdr,'Private_0019_100a');
      slc = double(hdr.Private_0019_100a);
    else
      slc = 1;
      
    end;
  else
    error('not a dicom image... need to specify the number of slices');
  end;
end;

if (nargin<3)
  f = ceil(sqrt(slc));
  m = m/f; n = n/f;
  f1 = f; f2 = f;
else
  f1 = m_in; m = m/f1;
  f2 = n_in; n = n/f2; 
end;

if ( exist('trefold') ~= 2 )
  addpath('~/matlab/mathlib/tensor');
end;

a2 = trefold( a.', [ n*f2 f1 m ], 1 );
a2 = trefold( a2, [ n m f1*f2 ], 1 );
a2 = permute( a2, [2 1 3]);


if (m ~= n),
  o = max([ m n ]);
  p = min([ m n ]);
  a3 = zeros([ o o f1*f2]);
  a3( (o-p)/2 + (1:p),:,:) = a2;
else
  a3 = a2;
end;

a3 = a3(:,:,1:slc);

