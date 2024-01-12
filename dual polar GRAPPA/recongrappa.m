function [F,I,mask] = recongrappa( varargin );
% [F,I,mask ] = recongrappa(Wk, Sk, Y, [opt args]);
% 
% opt args:
%  'v', {0,1}              - verbosity off/on
%  'kernel', string        - grappa kernel size, x by y. e.g. '2x5', '4x1'
%  'N', (xy)-by-L matrix   - grappa reconstruction parameters
%  'acs', array of indices - index locations to use for ACS lines
%  'dk', delta-k           - distance between phase encode lines used by
%                            the GRAPPA kernel.  default: 2
% 
% This is an implementation of     
% "Generalized autocalibrating partially parallel acquisitions"
%   by M. A. Griswold, P. M. Jakob, et. al..
%  Mag. Reson. Med.   47(6):1202-1210.  Jun 2002

% Copyright 2004-2007  Scott Hoge  (shoge at ieee dot org)

%
% uncombined images are generated for each coil in the array, by
% filling in k-space lines.
%
%
% GRAPPA: data from all coils used to patch data in one coil
%
%  coil 1     o * o * - * o *
%  coil 2     o * o * - * o *
%  coil 3     o * o * - * o *
%  coil 4     o * o * x * o *
%
%  *: acquired line, o: unacquired line, -: auto-callibration line
%  x: reconstructed line
%
%  S_j (k_y - m\Delta k_y ) = \sum_{l=1}^L \sum_{b=0}^{N_b-1}
%                                  n(j,b,l,m) S_l(k_y - b A \Delta k_y)
%
% in a variable density sampling scheme, the center of k-space is fully
% sampled and can be used to determine the linear fit parameters.
%
% A   - acceleration factor
% N_b - number of blocks used in the reconstruction
% 
% N_b = 1 ----> GRAPPA = VD-AUTO-SMASH

global prnt

Wk = varargin{1};
s  = varargin{2};
Y  = varargin{3};

acs = []; N = []; verbose=0; kernel = '2x3'; dk = 2;
prnt = 0;   % set to '1' to output some debugging information
Nbias = [];
chi = 1e-8;

vararg = varargin(4:end);
for cnt=1:2:length(vararg),
  if strcmp( vararg{cnt}, 'v' ),
    verbose = cell2mat(vararg(cnt+1));
  elseif strcmp( vararg{cnt}, 'N' ),
    N = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'Nbias' ),
    Nbias = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'kernel' ),
    kernel = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'dk' ),
    dk = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'chi' ),
    chi = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'acs' ),
    acs = vararg{cnt+1};
  elseif (strcmp( vararg{cnt}, 'report' ) | strcmp( vararg{cnt}, 'r' )),
    prnt = vararg{cnt+1};
  end;
end;

if isempty(acs), 
  d = diff(Y);
  acs = Y(find(d == 1));
  % acs = acs(2:length(acs));
end;

if length(size(Wk)) == 3, 
  [Wm,Wn,Wc] = size(Wk);
else,
  Wm = Wk(1); Wn = Wk(2); Wc = Wk(3);
end;

%% zero-pad the k-space data
Sk = zeros( Wm, Wn, Wc );
Sk(Y,:,:) = s;

% set the kernel spacings based on 'kernel' string
nky = str2double(kernel(1));               % should be even
if (nky==1),
  indx0 = 1;
else,
  indx0 = floor( [-(nky-1):2:(nky-1)] * dk/2 );
end;
% if (nky == '2'), indx0 = [-1 1];
% elseif (nky == '6'), indx0 = [ -5 -3 -1 1 3 5];
% else % (nky == '4'), indx0 = [ -3 -1 1 3];
% end;

nkx = str2double(kernel(3));               % should be odd
indx3 = ceil((0-nkx)/2) -1 + (1:nkx); % [-2:2];

rptupd( sprintf(['kernel size: ' num2str(length(indx0)) 'x' num2str(length(indx3)) ';']) );
rptupd( sprintf([' indx0 = [' sprintf(' %d',indx0) ']; ']) );
rptupd( sprintf([' indx3 = [' sprintf(' %d',indx3) '];\n']) );

stime = clock;

plist = zeros( length(acs), 1 );

for cnt2 = 1:length(acs);
  %% check if the GRAPPA pattern...
  indx = acs(cnt2) + indx0;
  
  %% can be accomodated by the sampling pattern 
  tmp = find( Y(:)*ones(1,length(indx)) == ones(length(Y),1)*indx );
  plist(cnt2) =  ( length(tmp) == length(indx) );
end;
%% use only those ACS lines that fit the GRAPPA kernel
acs = vec( acs(logical(plist)) ).';
rptupd([ ' ACS lines:' sprintf(' %d ',acs) '\n' ]);

rptupd('processing coil: ');
for l=1:size(s,3),
  rptupd([' ' num2str(l)]);
  %% determine the k-space filling coefficients for each coil:

  %%% setup the system matrix
  A2 = zeros( size(s,2), length(indx3)*length(indx0)*size(s,3) );

  if size(N,2) < l,                     % if coefficients are not supplied
  %%% Griswolds way:
  %   tic;
  %   A2 = zeros( size(s,2), length(indx3)*length(indx)*size(s,3) );% 20*size(s,3) );
  %   for x0=1:size(s,2),
  %     indx2 = x0 + indx3;
  %     indx2( indx2 < 1 ) = indx2( indx2 < 1 ) + size(s,2);
  %     indx2( indx2 > size(s,2) ) = indx2( indx2 > size(s,2) ) - size(s,2);
  %     
  %     A2(x0,:) = vec( permute(squeeze( Sk(indx,indx2,:)),[2 3 1 ]) ).';
  %   end;
  %   toc
  
  %%% My way: (heh, 8x faster)
  % tic;
  A = zeros( length(acs)*size(s,2), length(indx3)*length(indx0)*size(s,3) );
  A2 = zeros( size(s,2), length(indx3)*length(indx0)*size(s,3) );
  if (l==1), rptupd(['\n size(A) = ' num2str(size(A)) '\n']), end;

  for cnt2=1:length(acs),
    indx = acs(cnt2) + indx0;
    for cnt=1:length(indx3),
      y0 = mod( indx3(cnt) + [1:size(Sk,2)] , size(Sk,2) );
      y0( y0 == 0 ) = size(Sk,2);
      
      A( size(s,2)*(cnt2-1) + (1:size(s,2)), ...
         cnt:length(indx3):size(A,2) ) = ...
          reshape( permute(squeeze( Sk(indx,y0,:) ),[ 2 3 1 ]), ...
                   size(Sk,2), size(Sk,3)*length(indx) );
    end;
  end;
  % toc

  %%% setup the objective vector
  b = vec( Sk( acs, :, l ).' );
  
  %%% find the coefficients
  % n = A \ b; % pinv( A, norm(A)*0.01 ) * b;  
  
  if (size(Nbias,1)==0), % isempty('Nbias'),
    % n = cgsolv( A'*A, A'*b, zeros(size(A,2),1), size(A,2) );
    [ n tt ] = bicg( A'*A + chi*eye(size(A,2)), A'*b );
  else,
    n = cgsolv( A'*A, A'*b, Nbias(:,l), size(A,2)/4 );
  end;
  
  N(:,l) = n;                           % store them
  if verbose, pltcmplx( A*n, b ); keyboard; end;
  end;

  % extract the coefficients for this coil
  n = N(:,l);

  %% now go recon the missing lines of k-space:
  mss = zeros(Wm,1);            % set the missing lines index vector 
  mss(Y) = 1;

  for cnt0=1:length(mss),
    if (mss(cnt0) == 1);   continue; end;

    indx = cnt0 + indx0;
    if ( (min(indx) < 1) | (max(indx) > size(Sk,1))); continue; end;

    if ( sum(sum( Y(:)*ones(size(indx)) == ones(size(Y(:)))*indx )) < length(indx) ),
      continue; 
    end;

    %     for x0=1:size(s,2),
    %       indx2 = x0 ; % + [-2:2];
    %       indx2( indx2 < 1 ) = indx2( indx2 < 1 ) + size(s,2);
    %       indx2( indx2 > size(s,2) ) = indx2( indx2 > size(s,2) ) - size(s,2);
    %       
    %     end;
    for cnt=1:length(indx3),
      y0 = mod( indx3(cnt) + [1:size(Sk,2)] , size(Sk,2) );
      y0( y0 == 0 ) = size(Sk,2);
        
      A2( (1:size(s,2)), cnt:length(indx3):size(A2,2) ) = ...
          reshape( permute(squeeze( Sk(indx,y0,:) ),[ 2 3 1 ]), ...
                   size(Sk,2), size(Sk,3)*length(indx) );
    end;
    Sk( cnt0, :, l ) = A2 * n;
    % figure(1); imagesc( log10(abs(Sk(:,:,l)) + 1e-1) ); pause(0.5);
  end;

end;
rptupd('\n');
rptupd( sprintf( 'Elapsed time: %8.5g seconds.\n', etime(clock,stime) ) );

F = Sk;
I = sqrt( sum(abs(ifft2(Sk)).^2,3) ); 
mask = N;

clear global prnt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rptupd( str )
global prnt

if (prnt)
  fprintf( str );
end;

return;
