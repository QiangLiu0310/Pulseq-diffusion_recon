function [F,I,N] = recongrappa_multik( varargin );

% [F,I,N ] = recongrappa(Wk, Sk, Y, [opt args]);
% 
% opt args:
%  'v', {0,1}              - verbosity off/on
%  'kernel', string        - grappa kernel size, x by y. e.g. '2x5', '4x1'
%  'N', (xy)-by-num_patterns-L matrix   - grappa reconstruction parameters
%  'acs', array of indices - index locations to use for ACS lines
%  'radius', r             - reconstruction radius.  When GRAPPA is used to 
%                            improve data coverage for coil sensitivity
%                            estimates, only the center of k-space needs to
%                            be reconstructed. r \in [ 0 .. 1 ]
%  'iir', f                - option to control IIR-GRAPPA
%                            feedback. Default number of iir taps is 1/2
%                            the number of open lines in the kernel.  It
%                            can be set either higher, or set to zero (0)
%                            to turn IIR off.
% 
% This is an implementation of     
% "Generalized autocalibrating partially parallel acquisitions"
%   by M. A. Griswold, P. M. Jakob, et. al..
%  Mag. Reson. Med, 47(6):1202-1210, Jun 2002
%
% and expanded to use feedback taps, as in 
% "IIR GRAPPA for parallel MR image reconstruction."
%   by Chen,Zhang, Yang, Kellman, Johnston, and Egan.  
%   Magn Reson Med, 63(2):269-542, Feb 2010.

%
% Copyright 2004-2007  Scott Hoge  (shoge at ieee dot org)
%
% All rights reserved. Licensed according to the GNU General Public Licence
% (see http://www.gnu.org)
%

% uncombined images are generated for each coil in the array, by
% filling in k-space lines.
%
%
% AUTO-SMASH: estimate k-space lines via linear fit from neighboring lines
%
%  coil 1     o * o * a * o *
%  coil 2     o * o * a * o *
%  coil 3     o * o * a * o *
%  coil 4     o * o * a * o *
%
%  comp data  o * o * x * o *
%
% GRAPPA: data from all coils used to patch data in one coil
%
%  coil 1     o * o * a * o *
%  coil 2     o * o * a * o *
%  coil 3     o * o * a * o *
%  coil 4     o * o * x * o *
%
%  *: acquired line, o: unacquired line, a: auto-callibration line
%  x: reconstructed line
%
%  S_j (k_y - m\Delta k_y ) = \sum_{l=1}^L \sum_{b=0}^{N_b-1}
%                                  n(j,b,l,m) S_l(k_y - b A \Delta k_y)
%
% IIR taps have a simalar mathematical description, but can be drawn from
% both estimated and measured data. When IIR-GRAPPA is used, the recon
% starts from the ACS region and proceeds to the higher spatial frequency
% locations.  To maintain causality, two IIR pattern markers are used: 
% 'p' refers to IIR taps that have a higher index number than the target point,
% and 'm' refers to IIR taps that have a lower index number.
%
% in a variable density sampling scheme, the center of k-space is fully
% sampled and can be used to determine the linear fit parameters.
%
% A   - acceleration factor
% N_b - number of blocks used in the reconstruction
% 
% N_b = 1 ----> GRAPPA = VD-AUTO-SMASH

Wk = varargin{1};
s  = varargin{2};

if (nargin < 3)
  Y = [];
else
  Y  = varargin{3};
end;

if isempty(Y)
  Y = find( sum(abs(s(:,:,1)),2) );
  s = s(Y,:,:);
end;

acs = []; N = []; verbose=0; kernel = '2x3'; dks = []; prnt = 0;

iir = 0;
rdist = 1.0;

usenorm = 1;                            % use system normalization (eta)
eta = 1.0;
chi = 1e-6;
solver = 'bicg';

vararg = varargin(4:end);
for cnt=1:2:length(vararg),
  if strcmp( vararg{cnt}, 'v' ),
    verbose = cell2mat(vararg(cnt+1));
  elseif strcmp( vararg{cnt}, 'N' ),
    N = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'dks' ),
    dks = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'kernel' ),
    kernel = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'acs' ),
    acs = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'radius' ),
    rdist = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'normalize' ),
    usenorm = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'eta' ),
    eta = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'chi' ),
    chi = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'solver' ),
    solver = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'iir' ),
    iir = vararg{cnt+1};
  elseif (strcmp( vararg{cnt}, 'report' ) || strcmp( vararg{cnt}, 'r' )),
    prnt = vararg{cnt+1};
  end;
end;

if length(size(Wk)) == 3, 
  [Wm,Wn,ncoils] = size(Wk);
else,
  Wm = Wk(1); Wn = Wk(2); ncoils = Wk(3);
end;

if ~( (usenorm == 0 ) || ( usenorm == 1) )
  wrnstr = 'invalid ''normalization'' switch declared: ';
  if ischar(usenorm)
    wrnstr = [wrnstr '(' usenorm ')' ];
  else
    wrnstr = [wrnstr '(' num2str(usenorm) ')' ];
  end;
  warning([ wrnstr '. resetting to usenorm = 1' ]);
  usenorm = 1;
end;
  if isempty(chi)
    nlams = 20;
  else
    nlams = length(chi);
  end;

%%% correct for the toggling in acquisition...
% s = tmult( s, diag( (-1).^Y ), 1 );

if nargin < 3,
  Sk = s;
else
  %% zero-pad the k-space data
  Sk = zeros( Wm, Wn, ncoils );
  Sk(Y,:,:) = s;
end;

rptupd(prnt,['* ' mfilename ': kernel size = ' kernel '\n']);

d = diff(Y(:));
if ( isempty(acs) ),
  prfl = max( sum(abs( Sk(:,:,:)),3).' );

  if length( d == 1) == 1,
    acs = find( prfl == max(prfl(Y(find(d==1)))) );
  else
    acs = Y([ min(find(d==1)); find( d == 1)+1 ]);
  end;
  rptupd(prnt,['* ' mfilename ': acs =']);
  rptupd(prnt,sprintf(' %d', acs) );
  rptupd(prnt,'\n');
end;

if ( isempty(dks) ),
  %% determine which delta-k's to repair.
  dks = [];
  for cnt=max(d):-1:2,
    if (sum( cnt == d ) > 0),
      dks = [ cnt; dks ];
    end;
  end;
end;

% set the kernel spacings based on 'kernel' string
nky = str2double(kernel(1));           
nkx = str2double(kernel(3:end)); % should be odd

indx3 = ceil((0-nkx)/2) -1 + (1:nkx); % [-2:2];

pattern = {};

if ~isempty(N),
  % if N is supplied, set how many iir taps are used
  for cnt=1:length(N)
    pattern{cnt} = N(cnt).pattern;
  end;
  allptn = [ N(:).pattern ];
  nump = sum(allptn == 'p');
  numm = sum(allptn == 'm');
  if ( nump == numm ),
    iir = nump / length(N);
  else
    iir = 0
  end;
end;

if isempty(pattern)
%% foreach dk spacing
for cnt=1:length(dks),
  tmp = [ repmat([ '*' repmat('o',1,dks(cnt)-1) ],1,nky-1) '*' ];
  z2 = strfind(tmp,'o');
  if isempty(z2), 
    tmp = [ 'o' tmp 'o' ];   
    z2 = strfind(tmp,'o');
  end;
  [tmp2,indx] = sort( abs(z2 - z2( floor(median(1:length(z2))) )) );

  olength = length(pattern);
  for cnt2 = 1:length(z2), % (dks(cnt)-1), % 
    cnt3 = indx(cnt2);
    pattern{length(pattern)+1} = tmp;
    pattern{length(pattern)}(z2(cnt3)) = 'x';

    if isempty(iir),
      iirmxnum = round(sum(find( tmp == 'o' )~=0)/2);
      iir = iirmxnum;
    else, 
      iirmxnum = iir;
    end;
    if (iir~=0)
      %%% tag the leading-edge and following-edge locations, so that 
      %%% IIR-GRAPPA operates in a causal manner.

      iirndx = find( pattern{end} == 'o' );
      r = find( pattern{end} == 'x' );
      pattern{length(pattern)+1} = pattern{length(pattern)};

      pattern{end}( iirndx(iirndx > r) ) = 'p';
      if ( length(find( pattern{end} == 'p')) ~= iirmxnum )
        pndx = find( pattern{end} == 'p');
        if ( length(pndx) < iirmxnum )
          pattern{end}( end + [1:(iirmxnum-length(pndx))] ) = 'p';
        else
          pattern{end}( pndx ) = 'o';
          pattern{end}( pndx(1:iirmxnum) ) = 'p';
        end;
      end;
      
      pattern{end-1}( iirndx( iirndx < r ) ) = 'm';
      if ( length(find( pattern{end-1} == 'm')) ~= iirmxnum )
        pndx = find( pattern{end-1} == 'm');
        if ( length(pndx) < iirmxnum )
          pattern{end-1} = [ repmat('m',1,iirmxnum-length(pndx)) pattern{end-1} ];
        else
          pattern{end-1}( pndx ) = 'o';
          pattern{end-1}( pndx( end+1-(1:iirmxnum)) ) = 'm';
        end;
      end;
    end;
  
  end;

end;
end;

% pattern{1} = '*xoo*o*';
% pattern{2} = '*oxo*o*';
% pattern{3} = '*oox*o*';

% find the recon patterns we can use 
plist = cell(length(pattern),1);

for cnt = 1:length(pattern),
  v = extrapolate_pattern(pattern{cnt});     

  % if cnt == 14, keyboard; end;

  for cnt2 = 1:length(acs);
    %% check if the GRAPPA pattern...
    indx = acs(cnt2) + v;

    %% can be accomodated by the sampling pattern 
    tmp = find( Y(:)*ones(1,length(indx)) == ones(length(Y),1)*indx );
    if ( length(tmp) == length(indx) ),
      plist{cnt} =  [ plist{cnt} acs(cnt2) ];
    end;
  end;
end;

rptupd(prnt,['* ' mfilename ': sampling patterns in use = \n' ]);
for cnt1=1:length(pattern);
  rptupd(prnt,sprintf('  %s: %d acs lines =', pattern{cnt1}, length(plist{cnt1})) );
  rptupd(prnt,sprintf(' %3d', plist{cnt1}) );
  % for cnt2 = 1:length( plist{cnt1} ),
  %   rptupd(prnt,('  %s', pattern{plist{cnt1}} );
  %   % v = extrapolate_pattern(pattern{plist{cnt1}(cnt2)});
  %   % norm(vec(Sk( acs(cnt1)+v, :, : )));
  % end;
  rptupd(prnt,'\n');
end;
if (length(plist) == 0), 
  warning(['sorry, acquisition pattern is incompatible with available GRAPPA recon patterns']);
end;

st = clock; % tic;

if verbose, rptupd(prnt,'\n'); end;

mss = zeros( size(Sk,1), 1 );         % set the missing lines index vector 
mss(Y) = 1;

if (rdist>1.0) rdist = 1; end;
if (rdist<0.2) rdist = 0.2; end;
rd = floor(rdist*size(Sk,2));
kx_indx = (size(Sk,2)-rd)/2+(1:rd);

if isempty(N),
  % N = zeros( length(indx3)*nky*size(s,3),length(pattern),size(s,3));

  if exist('cgsolv') ~= 2,
    %% add a path to the conjugate gradient solver, if needed.
    ans = which('recongrappa_multik');
    idx = findstr('/',ans)-1;
    addpath([ ans(1:idx(end-2)) '/mathlib/conjgrad' ]);
  end;
  
  rptupd(prnt,'calc recon param\n');
  for cnt = 1:length(pattern), 
    v = extrapolate_pattern(pattern{cnt});
    f = identify_feedback_indices(pattern{cnt});
    nfc = length(f);
    
    %% and for each pattern ... 
    %% A2 = [];     b2 = [];
    
    %%  A is data from the expanded input data matrix.
    %%     size(A) = [ kxnum * nacs, ncoils * nky * nkx ]
    %%
    %% in the matlab (and c-code) structure is:
    %%    for each element of indx3
    %%                   coil 1   |    coil 2   |    coil 3   | ...
    %%                -3 -1  1  3 | -3 -1  1  3 | -3 -1  1  3 | ...
    %%                 o  o  o  o |  o  o  o  o |  o  o  o  o | ...
    %%             ^
    %%  ACS(1) xres|
    %%             v
    %%             ^
    %%  ACS(2) xres|
    %%             v
    %%             ^
    %%  ACS(3) xres|
    %%             v
    %% 
    %%  the older version of matlab code sorted the columns by
    %% 
    %%     | c1 c2 c3 ... | c1 c2 c3 ... | c1 c2 c3 ...  | c1 c2 c3 ...
    %%     | -3 -3 -3 ... | -1 -1 -1 ... |  1  1  1 ...  |  3  3  3 ...
    %%     |  o  o  o ... |  o  o  o ... |  o  o  o ...  |  o  o  o ...
    %% 
    %%  the order is immaterial, as long at it is consistent...
    %%  UPDATE (12 Jan 2011): the order does matter for GROG, so we have
    %%                        now changed it here so that the c-code and matlab code match.
    %%   
    %%  The old way to switch between the two was:
    %%    Accode(:, permute( reshape(1:size(Accode,2),YKS,Ncoils,XKS), [3 2 1]) ) == Amatlab
    %%    Accode = Amatlab(:, permute( reshape(1:size(Accode,2),XKS,Ncoils,YKS), [3 2 1]))
    A = zeros( length(plist{cnt})*length(kx_indx), length(indx3)*nky*size(s,3) );
    A_bak = zeros( length(plist{cnt})*length(kx_indx), length(indx3)*nky*size(s,3) );
    B = zeros( length(plist{cnt})*length(kx_indx), length(indx3)*nfc*size(s,3) );
    
    rptupd(prnt,sprintf('%3d  %s: ', cnt, pattern{cnt}));

    % v = extrapolate_pattern(pattern{cnt});
    N(length(N)+1).pattern = pattern{cnt};

    n = zeros( (nky+nfc)*size(s,3)*length(indx3),size(s,3),nlams);

    p_src = ones([1 length(plist{cnt}) ]);
    for cnt2 = 1:length(plist{cnt}),
      %% and ACS line ...

      %% if the sampling pattern accomodates the GRAPPA pattern...
      indx  = plist{cnt}(cnt2) + v;         % feedforward indices
      indx2 = plist{cnt}(cnt2) + f;         % feedback indices

      %% determine the k-space filling coefficients 
      for cnt3=1:length(indx3),
        y0 = mod( indx3(cnt3) + [kx_indx] , size(Sk,2) );
        y0( y0 == 0 ) = size(Sk,2);
        
        % A( size(s,2)*(cnt2-1) + (1:size(s,2)), cnt3:length(indx3):size(A,2) ) = ...
        A( rd*(cnt2-1) + (1:rd), (cnt3-1)*nky*ncoils + (1:nky*ncoils) ) = ...
            reshape( permute( Sk(indx,y0,:), [ 2 1 3 ]), rd, size(Sk,3)*length(indx) );

        % B( size(s,2)*(cnt2-1) + (1:size(s,2)), cnt3:length(indx3):size(B,2) ) = ...
        B( size(s,2)*(cnt2-1) + (1:size(s,2)), (cnt3-1)*length(indx2)*ncoils + (1:length(indx2)*ncoils) ) = ...
            reshape( permute( Sk(indx2,y0,:), [ 2 1 3 ]), rd, size(Sk,3)*length(indx2) );
      
      end;
      %% for an N-by-1 kernel, these next lines are enough...
      % A2 = [ A2; ...
      %         reshape( permute(squeeze( Sk(indx,:,:) ),[ 2 3 1 ]), ...
      %                  size(s,2), size(s,3)*length(indx) )
      %       ];
      % 
      % b2 = [b2; squeeze( Sk( plist{cnt}(cnt2), :, l ) ).' ];
    end;

    % normalize, according to Algo 4.2 of M. Schneider, "GPGPU for Accelerated GRAPPA Autocalibration in
    %    Magnetic Resonance Imaging," Masters Thesis, U of Erlangen, Apr 2008
    % A_bak = A;
    if (usenorm)
      p_src = sum(abs(A),2).^(-eta/2);    % compute the power in each line of the linsys
      % savend([ 'Araw_g_' pattern{cnt} '.nd' ], A);
      Abak = A;
      A = tmult(A,diag(p_src),1);         % normalize by the power in each linsys
    end;
    % Apinv = pinv(A);
    % AA = A'*A;
    AA = [A B]'*[A B];
    % if (usenorm)
    %   AA = AA + diag(norm(AA)*chi);
    % end;
    At = [A B]';

    btmp = Sk( plist{cnt}, :, : );
    for l=1:size(Sk,3),
      rptupd(prnt,'.');
      % A is the same for each coil ...
      b = vec( btmp(:,:,l).' );

      if (usenorm)
        b = b.*p_src(:);                  % normalize by the power in each linsys
      end;

      % n(:,l) = A \ b;
      % n(:,l) = Apinv * b;
      % n(:,l) = cgsolv( AA, At*b, zeros(size(At,1),1), size(At,1) );
      if length(chi) == 1,
        if strcmp(solver,'cgsr')
          n(:,l) = cgsolv( AA + chi*eye(size(AA)), At*b, zeros(size(At,1),1), size(At,1) );
        elseif strcmp(solver,'lsqr')
          [n(:,l),~,~] = lsqr_hybrid(AA, At*b, chi, 20);
        else
          [n(:,l), ~ ] = bicg( AA + chi*eye(size(AA)), At*b );
        end;
      else
        [n(:,l,:),~,~] = lsqr_hybrid(AA, At*b, chi, 20);
      end;
      
      if (verbose>2), pltcmplx( tmult( n(:,l,:), At', 1), repmat(b,[1 1 length(chi)]) ); keyboard; end;
    end;
    N(length(N)).n = n;
    N(length(N)).chi = chi;
    N(length(N)).eta = eta;

    rptupd(prnt,'\n');
  end;  
end;

rptupd(prnt,'recon missing lines\n');
%% now go recon the missing lines of k-space:

if (verbose); keyboard; end;

kypts = 1:length(mss);
tmp = round(mean(kypts)-1);
kypts( 1:tmp ) = kypts( tmp:-1:1 );

linerpt = repmat(' ',1,length(mss));

for kycnt=1:length(mss),
  cnt2 = kypts(kycnt);
  
  if verbose, rptupd(prnt,sprintf(' line %3d, from lines ',cnt2)); end;

  if ( find( Y == cnt2) ); 
    % rptupd(prnt,'*'); 
    linerpt(cnt2) = '*';
    if verbose, rptupd(prnt,'\n'); end;
    continue; 
  end;

  tmp = [];
  %% determine which recon pattern we can use:
  for cnt4 = 1:length(pattern),
    v = extrapolate_pattern(pattern{cnt4});
    f = identify_feedback_indices(pattern{cnt4});
    nfc = length(f);

    if ( (kycnt < length(mss)/2) && sum(find( pattern{cnt4} == 'm')) ), continue, end;
    if ( (kycnt > length(mss)/2) && sum(find( pattern{cnt4} == 'p')) ), continue, end;

    indx  = cnt2 + v;
    indx2 = cnt2 + f;
    %     tmp = find( Y(:)*ones(1,length(indx)) == ones(length(Y),1)*indx );
    %     if  ( ... % ( (min(indx) > 1) & (max(indx) < size(Sk,1)) ) & ...
    %         ( length(tmp) == length(indx) ) ),
    %       break;                          % jump out of the loop
    %     end;
    if ( indx(1) < 1) | (indx(end) > size(Sk,1)); continue; end;
    
    tmp = (mss( indx ) == ones(length(indx),1));
    if ( sum(tmp) == length(tmp)  ),
      break; 
    end;

  end;  
  
  if ( ( (min(indx) < 1) | (max(indx) > size(Sk,1)) ) | ...
       ( sum(tmp) ~= length(indx) ) ),
    % 
    linerpt(cnt2) = 'o';
    if verbose, rptupd(prnt,'o'); rptupd(prnt,'\n'); end;
    continue; 
  end;
  
  %% recon k-space line based on that pattern:
  A = zeros( size(s,2), length(indx3)*length(v)*size(s,3) );
  B = zeros( size(s,2), length(indx3)*length(f)*size(s,3) );
  
  for cnt3=1:length(indx3),
    y0 = mod( indx3(cnt3) + [1:size(Sk,2)] , size(Sk,2) );
    y0( y0 == 0 ) = size(Sk,2);
    
    % A( (1:size(s,2)), cnt3:length(indx3):size(A,2) ) = ...
    A( (1:size(s,2)), (cnt3-1)*nky*ncoils + (1:nky*ncoils) ) = ...
        reshape( permute( Sk(indx,y0,:), [ 2 1 3 ]), ...
                 size(Sk,2), size(Sk,3)*length(indx) );

    % B( (1:size(s,2)), cnt3:length(indx3):size(B,2) ) = ...
    B( (1:size(s,2)), (cnt3-1)*length(indx2)*ncoils + (1:length(indx2)*ncoils) ) = ...
        reshape( permute( Sk(indx2,y0,:), [ 2 1 3 ]), ...
                 size(Sk,2), size(Sk,3)*length(indx2) );
  end;

  % A2 = reshape( permute(squeeze( Sk(indx,:,:) ),[ 2 3 1 ]), ...
  %              size(Sk,2), size(Sk,3)*length(indx) );

  if (verbose), % ( sum(Sk(Y(mss(cnt2))+1,:,l)) ~= 0 ); 
    rptupd(prnt,sprintf(' %3d',[ indx(:) ]));
    rptupd(prnt,sprintf(' (pattern: %s)', pattern{cnt4}) );
    rptupd(prnt,'\n');
    if (verbose>1),
      figure(1); imagesc( log10(abs(Sk(:,:,1)) + 1e-1) ); pause(0.5);
      keyboard;
    end;
  end;
  linerpt(cnt2) = 'x';

  for cnt3=1:length(N),
    if strcmp( N(cnt3).pattern, pattern{cnt4} ),
      break;                            % jump out of the loop
    end;
  end;
  if isfield( N(cnt3),'coef')
    n = N(cnt3).coef;
  else
    n = N(cnt3).n;
  end;
  
  for l=1:size(Sk,3),
    % for each coil ...

    % if cnt2==109, keyboard; end;
    % rptupd(prnt,(' line %3d, pattern %3d\n', cnt2, plist(cnt3) );
    Sk( cnt2, :, l ) = [A B] * n(:,l);
  end;

  if (verbose>1), keyboard; end;
end;
rptupd(prnt,linerpt);
rptupd(prnt,'  \n');

% Sk(1:2:size(Sk,1),:,:) = -Sk(1:2:size(Sk,1),:,:);
% Sk(:,1:2:size(Sk,2),:) = -Sk(:,1:2:size(Sk,2),:);

F = Sk;
I = sqrt( sum(abs(ifft2(Sk)).^2,3) ); 
I = I./(prod(size(I)));
mask = 0;

rpttime = etime(clock,st);
rptupd(prnt,sprintf(' GRAPPA recon time: %g (avg per kernel: %g)\n',  ... 
        rpttime, rpttime/length(pattern)) ); 

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rptupd( prnt, str )

if (prnt)
  fprintf( str );
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = extrapolate_pattern(p)
 
r = find( p == 'x' );
v = find( p == '*' ) - r;
return;

function f = identify_feedback_indices(p)
 
r = find( p == 'x' );
tmp = (p == 'p') + (p == 'm') ;
if (sum(tmp) ~= 0),
  f = find( tmp ) - r;
else,
  f = [];
end;

return;

