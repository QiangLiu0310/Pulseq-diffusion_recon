function varargout = phzshift(a,b,mode,verbose,phz,thrsh);
% [c_out,y_out] = phzshift( a, b, [mode, verbose, phz]);
%
% calculates the coordinate shift between two k-space data sets, a and b.
% Using the Phase Correlation Method, each k-space set is transformed
% to the image domain, the phase difference is calculated and then
% applied to the second set.  After this correction, the two data sets
% are added together and output.
%
% inputs can be any number of dimensions, but 'a' and 'b' must be the
% same size
%
% c_out is the combination of (a+b)/2, after phase correction.
% 
% y_out contains the estimated phase shift between 'a' and 'b' for
% each frame.
%
% Advanced (optional) options:
%    mode - a cell variable containing flags to various settings.
%           'nofft'    : The input data is already in the spatial domain
%           'nocombo'  : Returns  b with the phase correction applied.
%           'nofilter' : use the full data to estimate the phase shift, 
%                        instead of the default setting which uses a
%                        low-passed filtered version.
%           'center'   : center the low-pass filter.
% verbose - set to '1' to show debugging information
%     phz~ - if nonempty, use this matrix instead of estimating the phase
%           shift.  (i.e. y_out from a previous run)
%   thrsh - the mask threshold for the svd calculation. range: [0:1]

% 2013/nov/07: jonathan polimeni <jonp@nmr.mgh.harvard.edu>
% fixed non-integer index, fixed missing conditional for "nofilter" option

% 2014/sep/25: jonathan polimeni <jonp@nmr.mgh.harvard.edu>
% replaced phase difference fitting with Local Phase Correction (LPC) method
% based on the methods described in:
%
%   Feiweier T. Magnetic Resonance Method and Apparatus to Determine Phase
%   Correction Parameters. U.S. Patent No. 8,497,681, Issued 7/30/2013.

  
[m,n,c]=size(a);
c_out = zeros(m,n,c);
if nargout > 1,
  y_out = zeros(2,c); % coil channels why double the size?
end;

if nargin < 3,
  mode = {''};
end;

if nargin < 4,
  verbose = 0;
end;

if nargin < 5,
  phz = [];
end;
if ~isempty(phz),
  phz = reshape(phz,[2 c]);
end;

if nargin < 6,
  thrsh = 0.2;
end;


if ( sum(strcmp(mode,'nofilter')) || sum(strcmp(mode,'nofft')) ),
  f1 = 1; f2 = 1;
else
  % define a low-pass envelope to apply in kspace
  
  alph = 20;
  
  if (alph > m), alph2 = m; else, alph2 = alph; end;
  f1 = [ zeros(floor((m-alph2)/2),1); gausswin(alph2); zeros(ceil((m-alph2)/2),1)]; 

  if (alph > n), alph2 = n; else, alph2 = alph; end;
  f2 = [ zeros(floor((n-alph2)/2),1); gausswin(alph2); zeros(ceil((n-alph2)/2),1)]; 

  [tmp,maxx] = max(max(sum(abs(a(:,:,:)),3) ));
  [tmp,maxy] = max(max(sum(abs(a(:,:,:)),3)'));

  %%% correct the window coordinates to set ontop of DC (which may be off center)
  [ans,j2]=max(f1);
  [ans,j3]=max(f2);
  
  if sum( strcmp(mode,'center') )
    c1 = 0; c2 = 0;
  else,
    c1 = j2 - maxy;
    c2 = j3 - maxx;
  end;
  % fprintf('c1: %d c2: %d\n',c1,c2);
  
  if c1>0;
    f1 = f1( [ c1:length(f1) 1:(c1-1)]);
  elseif c1<0
    f1 = f1( [ length(f1)+(c1:0) 1:(length(f1)+c1-1)]);
  end;
  
  if c2>0;
    f2 = f2( [ c2:length(f2) 1:(c2-1)]);
  elseif c2<0
    f2 = f2( [ length(f2)+(c2:0) 1:(length(f2)+c2-1)]);
  end;

end;


% convert k-space to image-space
if sum(strcmp(mode,'nofft'));
  A = a;
  B = b;
elseif sum(strcmp(mode,'nofilter')),
  A = fftshift(ifft2(fftshift( reshape(a,[m n c]) )));
  B = fftshift(ifft2(fftshift( reshape(b,[m n c]) )));
else,
  A = fftshift(ifft2(fftshift( repmat( (f1*f2'),[1 1 c]).*reshape(a,[m n c]) )));
  B = fftshift(ifft2(fftshift( repmat( (f1*f2'),[1 1 c]).*reshape(b,[m n c]) )));
end;

mask = ( abs(A) > thrsh*max( abs(A(:)) ) );

A1 = reshape( permute(A,[2 1 3 4]),[n m*c]).';  % FE PE*COIL
B1 = reshape( permute(B,[2 1 3 4]),[n m*c]).';
M1 = reshape( permute(mask,[2 1 3 4]),[n m*c]).';

if (1),
  if isempty(phz)
    % calculate the normalized cross-correlation
    C = zeros(size(A1));
    % [JRP] why does next definition of "tt" below use a threshold of 0.1*max?
    tt = find( abs(A1(:)) > 0.01*max(abs(A1(:))) );
    C(tt) = A1(tt).*conj(B1(tt))./abs(A1(tt).*conj(A1(tt)));
  
    % extract out the dominant x and y shift vectors
    [u d v]=svds(M1.*C,1);
    % v = rsvd( mask.*C, 1 );
  

    m2 = sum(M1,1).';
    t2 = find(m2~=0);

    % [JRP] previously the vector "v" was unwrapped after extracting a
    % subset defined by vector "t2", but since the subset of samples may be
    % non-contiguous the unwrapping could add erroroneous phase
    % jumps---better to unwrap the vector first then extract the subset
    % afterwards
    unv = unwrap(angle(v));
    
    % calculate the shift from the slope of the phase.
    % here, we are only concerned with shifts along rows
    if (verbose),
      fprintf(' %d ',t2);
      fprintf('\n');
      Atmp = [ ones(length(t2),1) t2(:)-size(B,2) ];
      fprintf('A''*A  =');
      fprintf(' %d ', Atmp'*Atmp );
      fprintf('\n');
      fprintf('A''*v2 =');
      fprintf(' %f ', Atmp'*unv(t2) );
      fprintf('\n');
    end;
    if (0),
      % [JRP] "robustfit" performs better than backslash when outliers are
      % present and are not properly masked
      %y2 = ( [ ones(length(t2),1) t2(:)-size(B,2) ] ) \ unv(t2);
      y2 = robustfit(t2(:)-size(B,2), unv(t2) );
    else,
      % [JRP] based on "Local Phase Correction" of Feiweier
      dphz = angle(sum(v(2:end) .* conj(v(1:end-1))));
      y2 = [1; dphz];
    end;
    if (verbose) % verbose
      figure; plot( t2-size(B,2), unv(t2), 'bx-', t2-size(B,2), [ ones(length(t2),1) t2(:) ]*y2,'k-.' );

      fprintf(' x = %f  %f \n', y2 );
      keyboard;
    end;
  else,
    y2(2) = phz(1);
  end; % if exist('phz')
  
  % correct the phase of the second data set to match the first
  % (i.e. shift the k-space grids to match)
  if sum(strcmp(mode,'nofft'));
    A = a;% (:,:,cnt);
    B = b;% (:,:,cnt);
  else,
    A = ffts(ifft2(ffts( a,1:2 )), 1:2 );
    B = ffts(ifft2(ffts( b,1:2 )), 1:2 );
  end;
  for cnt=1:c,
    if sum(strcmp(mode,'nocombo'))
      B(:,:,cnt) = B(:,:,cnt) * diag( exp( -j*y2(2)*( (1:size(B,2))-size(B,2)/2 ) ) );
    else 
      A(:,:,cnt) = A(:,:,cnt) * diag( exp( +j*y2(2)/2*( (1:size(A,2))-size(A,2)/2 ) ) );
      B(:,:,cnt) = B(:,:,cnt) * diag( exp( -j*y2(2)/2*( (1:size(B,2))-size(B,2)/2 ) ) );
   end;
  end;
%  m1 = sum(mask,2).';
%  t1 = find(m1~=0);
%  % calculate the shift from the slope of the phase.
%  % here, we are only concerned with shifts along rows
%  y1 = ( [ ones(length(t1),1) t1(:) ] ) \ unv(t1);
%  
%  % correct the phase of the second data set to match the first
%  % (i.e. shift the k-space grids to match)
%  B = diag( exp( j*y1(2)*((1:size(B,1))-1) ) ) * B;

  % calculate the scalar phase difference between the sets ...
  if isempty(phz)
    if (0),
    P = zeros(size(A));
    % [JRP] why does previous definition of "tt" above use a threshold of 0.01*max?
    tt = find(A(:)>0.1*max(A(:)));
    P(tt) = A(tt).*conj(B(tt))./abs(A(tt).*conj(A(tt)));
    indx1 = round((size(P,1) - 10)/2) + (1:10); if ( sum(indx1<=0) > 0 ) indx1 = 1:size(P,1); end;
    x = mean(vec(P(  indx1, (end-10)/2 +(1:10), : ) ));
    else,
      % [JRP] based on "Local Phase Correction" of Feiweier
      x = sum(vec(A.*conj(B))); 
    end;
  else,
    x = phz(2);
  end;

  % and apply a scalar correction as well.
  B = B.*exp(j*angle(x));

  if nargout > 1,
    y_out(:,cnt) = [ y2(2); exp(j*angle(x)) ];
    if (verbose)
      fprintf(' p = %f  %f + j %f\n',y_out(1,cnt), ...
      real(y_out(2,cnt)),imag(y_out(2,cnt)) );
    end;
  end;

  % add the two sets together, and convert back to k-space for output.
  if sum(strcmp(mode,'nocombo'))
    c_out = B;
  else,
    c_out = (A+B)/2;
  end;
  if sum(strcmp(mode,'nofft')),
  else,
    c_out = iffts(fft2(iffts(c_out,1:2)),1:2);
  end;

end;


c_out = reshape( c_out, size(a) );
sz = size(c_out);
if ( (nargout > 1) && length(sz)>2),
  y_out = reshape( y_out, [ 2 sz(3:end) ] );
end;

% return output
if nargout > 0
  varargout{1}= c_out;
end

if nargout > 1
  varargout{2}= y_out;
end
 
if nargout > 2
  varargout{3}= y2;
end

if nargout > 3
  varargout{4}= x;
end
