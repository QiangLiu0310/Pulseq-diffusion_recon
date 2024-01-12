function [F,I,N] = dpg_sms( varargin );

% [F,I,N ] = dpg_segepi(Wk, Sk, Y, [opt args]);
% 
% opt args:
%  'v', {0,1}              - verbosity off/on
%  'kernel', string        - grappa kernel size, x by y.
%                            e.g. '2x5', '4x1', '3x5', '1x5'
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
% in a variable density sampling scheme, the center of k-space is fully
% sampled and can be used to determine the linear fit parameters.
%
% A   - acceleration factor
% N_b - number of blocks used in the reconstruction
% 
% N_b = 1 ----> GRAPPA = VD-AUTO-SMASH

Wk  = varargin{1};
kin = varargin{2};
Y   = varargin{3};

acs = []; N = {}; verbose=0; kernel = '2x5'; dks = []; prnt = 0;

iir = 0;
rdist = 1.0;
solver = 'bicg';
fast = 0;
seg = 1;

usenorm = 1;                             % use system normalization (eta)
eta = 1.0;
chi = 1e-6;

vararg = varargin(4:end);
for cnt=1:2:length(vararg),
  if strcmp( vararg{cnt}, 'v' ),
    verbose = cell2mat(vararg(cnt+1));
  elseif strcmp( vararg{cnt}, 'N' ),
    N = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'dks' ),
    dks = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'seg' ),
    seg = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'kernel' ),
    kernel = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'acs' ),
    acs = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'radius' ),
    rdist = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'normalize' ),
    usenorm = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'solver' ),
    solver = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'eta' ),
    eta = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'chi' ),
    chi = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'iir' ),
    iir = vararg{cnt+1};
  elseif (strcmp( vararg{cnt}, 'report' ) || strcmp( vararg{cnt}, 'r' )),
    prnt = vararg{cnt+1};
  elseif (strcmp( vararg{cnt}, 'fast' )),
    fast = 1;
  end;
end;

if length(size(Wk)) ==  4, 
  [Wm,Wn,Nslc,ncoils] = size(Wk);
else
  Wm = Wk(1); Wn = Wk(2); Nslc = Wk(3); ncoils = Wk(4);
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
% fprintf(' lin sys normalization? %d  (eta: %d, chi: %d)\n', usenorm,eta,chi);

%%% correct for the toggling in acquisition...
% s = tmult( s, diag( (-1).^Y ), 1 );

%% zero-pad the k-space data
dualks = 1;
Sk = zeros( Wm, Wn, Nslc, ncoils );
if iscell(kin)
  ROp = kin{1};
  ROn = kin{2};
  s = kin{3};
  clear kin;
elseif isstruct(kin);
  ROp = kin.p;
  ROn = kin.n;
  s = kin.target;
  clear kin;
else
  dualks = 0;
  s = kin; clear kin;
end;
Sk(Y,:,:,:) = s;

if ( sum( size(ROp) == size(Sk) ) ~= 4 )
  error('source and target data need to fill arrays of equal size');
end;

rptupd(prnt,['* ' mfilename ': kernel size = ' kernel '\n']);

d = diff(Y(:));
if ( isempty(acs) ),
  prfl = max( sum(sum(abs( Sk(:,:,:,:)),4),3).' );

  if length( d == 1) == 1,
    acs = find( prfl == max(prfl(Y(d==1))) );
  else
    acs = Y([ find(d==1, 1 ); find( d == 1)+1 ]);
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
  
  if ( dks(cnt) == 1 )
    % for no acceleration, just specify the pattern
    tmppat = repmat('+-',[1 nky]);
    tmppat = tmppat(1:nky);
    tmppat = vec(reshape( tmppat,[ nky/seg seg])')';
    for patcnt = 1:nky,
      pattern{length(pattern)+1} = tmppat;
      if ( pattern{length(pattern)}(patcnt) == '+' )
        pattern{length(pattern)}(patcnt) = '>';
      else
        pattern{length(pattern)}(patcnt) = '<';
      end;
    end;

    tmppat = repmat('-+',[1 nky]);
    tmppat = tmppat(1:nky);
    tmppat = vec(reshape( tmppat,[ nky/seg seg])')';
    for patcnt = 1:nky,
      pattern{length(pattern)+1} = tmppat;
      if ( pattern{length(pattern)}(patcnt) == '+' )
        pattern{length(pattern)}(patcnt) = '>';
      else
        pattern{length(pattern)}(patcnt) = '<';
      end;
    end;

    continue;
  end;

  % tmp = [ repmat([ '*' repmat('o',1,dks(cnt)-1) ],1,nky-1) '*' ];
  % tmp = repmat('o',1,dks(cnt)*seg+1);
  tmp = [];
  for nnn=1:nky/seg,
    tmp = [ tmp ...
            repmat([ '+' repmat('o',1,dks(cnt)-1) ],1,seg) ...
            repmat([ '-' repmat('o',1,dks(cnt)-1) ],1,seg) ];
  end;
  tmp = tmp( 1 : (dks(cnt)*nky - 1) );
  % tmp = [ tmp 'o' ];

  if ( nky/2 == floor(nky/2) ),
    tmp = tmp( 1 : find(tmp=='-', 1, 'last' ) );
  else
    tmp = tmp( 1 : find(tmp=='+', 1, 'last' ) );
  end;

  olength = length(pattern);
  if (dualks == 1)
    % for dual-kernels, we want to genereate coefficients for the sampled
    % lines as well.
    npat = 2*length(tmp);
  else
    npat = length(z2);
  end
  for cnt2 = 1:npat, % (dks(cnt)-1), % 

    if (dualks == 0 )
      cnt3 = indx(cnt2);
      pattern{length(pattern)+1} = tmp;
      pattern{length(pattern)}(z2(cnt3)) = 'x';
    else
      cnt3 = mod( cnt2-1, length(tmp) ) + 1;
      pattern{length(pattern)+1} = tmp;
      if (cnt2 > cnt3)
        % for the second round of patterns, swap the +/- locations
        indxp = find( tmp == '+' );
        indxm = find( tmp == '-' );
        pattern{length(pattern)}(indxp) = '-';
        pattern{length(pattern)}(indxm) = '+';
      end;

      if ( pattern{length(pattern)}( cnt3 ) == '+' )
        pattern{length(pattern)}( cnt3 ) = '>';
      elseif ( pattern{length(pattern)}( cnt3 ) == '-' )
        pattern{length(pattern)}( cnt3 ) = '<';
      else
        pattern{length(pattern)}( cnt3 ) = 'x';
      end;
%       if ( cnt2 <= npat/2 )
%         pattern{length(pattern)}( cnt3 ) = '<';
%       else
%         pattern{length(pattern)}( cnt3 ) = '>';
%       end;
    end;

  end;
end;
end;

% pattern{1} = '*xoo*o*';
% pattern{2} = '*oxo*o*';
% pattern{3} = '*oox*o*';

if (fast), 
  for cnt=1:length(pattern)/2;
    pattern2{cnt} = pattern{cnt};
  end;
  pattern = pattern2;
  clear pattern2;
end;

% find the recon patterns we can use 
plist = cell(length(pattern),1);

for cnt = 1:length(pattern),
  v = extrapolate_pattern(pattern{cnt},nky);     

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
  %   % v = extrapolate_pattern(pattern{plist{cnt1}(cnt2)},nky);
  %   % norm(vec(Sk( acs(cnt1)+v, :, : )));
  % end;
  rptupd(prnt,'\n');
end;
if (length(plist) == 0), 
  error(['sorry, acquisition pattern is incompatible with available GRAPPA recon patterns']);
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

  %   if exist('cgsolv') ~= 2,
  if exist('bicg') ~= 2,
    %% add a path to the conjugate gradient solver, if needed.
    ans = which('recongrappa_multik');
    idx = findstr('mri/parallel',ans)-1;
    % addpath([ ans(1:idx) 'mathlib/conjgrad' ]);
    addpath([ ans(1:idx) 'mathlib/templates' ]);
  end;
  
  rptupd(prnt,'calc recon param\n');
  for cnt = 1:length(pattern), 
    fprintf('-');
    v = extrapolate_pattern(pattern{cnt},nky);
%    f = identify_feedback_indices(pattern{cnt});
%    nfc = length(f);
    nfc = 0;

    p = find( (pattern{cnt}=='+') | (pattern{cnt}=='-') | (pattern{cnt}=='<') | (pattern{cnt}=='>') );
    if ( length(p) > nky )
      p = find( (pattern{cnt} == '+') | (pattern{cnt} == '-') );
    end;

    %indxpos = find( (pattern{cnt}(p) == '+') );
    %if ( length(indxpos) ~= nky/2 )
      % indxpos = find( (pattern{cnt}(p) == '+') | (pattern{cnt}(p) == '>') | (pattern{cnt}(p) == '<') );
      indxpos = find( (pattern{cnt}(p) == '+') | (pattern{cnt}(p) == '>') );
    %end;

    %indxneg = find( (pattern{cnt}(p) == '-') );
    %if ( length(indxneg) ~= nky/2 )
      % indxneg = find( (pattern{cnt}(p) == '-') | (pattern{cnt}(p) == '>') | (pattern{cnt}(p) == '<') );
      indxneg = find( (pattern{cnt}(p) == '-') | (pattern{cnt}(p) == '<') );
    %end;

    % fprintf('pattern: %s  ip:',pattern{cnt});
    % fprintf(' %d',indxpos); fprintf(' in:');
    % fprintf(' %d',indxneg); fprintf('\n');

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
    % A = zeros( length(plist{cnt})*size(s,2), length(indx3)*nky*size(s,3) );
    % B = zeros( length(plist{cnt})*size(s,2), length(indx3)*nfc*size(s,3) );
    for cnt4=1:Nslc,
      A{cnt4} = zeros( length(plist{cnt})*length(kx_indx), length(indx3)*nky*size(s,3) );

      % v = extrapolate_pattern(pattern{cnt});
      if (length(N) < cnt4); N{cnt4} = [] ; end;
      N{cnt4}(length(N{cnt4})+1).pattern = pattern{cnt};
    end;
    rptupd(prnt,sprintf('%3d  %s: ', cnt, pattern{cnt}));

    n = zeros( (nky+nfc)*ncoils*length(indx3),ncoils);

    for cnt2 = 1:length(plist{cnt}),
      %% and ACS line ...

      %% if the sampling pattern accomodates the GRAPPA pattern...
      indx  = plist{cnt}(cnt2) + v;         % feedforward indices
      % indx2 = plist{cnt}(cnt2) + f;         % feedback indices
      %% determine the k-space filling coefficients 
      for cnt3=1:length(indx3),
        % y0 = mod( indx3(cnt3) + [1:size(Sk,2)] , size(Sk,2) );
        y0 = mod( indx3(cnt3) + [kx_indx] , size(Sk,2) );
        y0( y0 == 0 ) = size(Sk,2);
        
        % generate a hybrid source matrix
        if (dualks)
          tmpSksub = zeros( size( Sk(indx,y0,:,:) ) );

          tmpSksub( indxpos, :, : ) = ROp( indx(indxpos),y0,: );
          tmpSksub( indxneg, :, : ) = ROn( indx(indxneg),y0,: );
        else
          tmpSksub = Sk(indx,y0,:,:);
        end;
        % A( size(s,2)*(cnt2-1) + (1:size(s,2)), cnt3:length(indx3):size(A,2) ) = ...
        for cnt4=1:Nslc,
          % build the system matrix for each slice
          A{cnt4}( rd*(cnt2-1) + (1:rd), (cnt3-1)*nky*ncoils + (1:nky*ncoils) ) = ...
              reshape( permute( (tmpSksub(:,:,cnt4,:)), [ 2 1 4 3 ]), rd, ncoils*length(indx) );
        end;
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
      for cnt4=1:Nslc,
        p_src{cnt4} = sum(abs(A{cnt4}).^2,2).^(-eta/2);    % compute the power in each line of the linsys
      end;
      %      if (cnt==7), keyboard; end;
    end;

    % Apinv = pinv(A);
    % AA = A'*A;
    AA = 0; 
    for cnt4=1:Nslc,                    % collapse the data matrix for the input slices
      if (usenorm)
        AA = AA + [ A{cnt4} ]'* diag(p_src{cnt4}.^2) *[ A{cnt4} ];
      else
        AA = AA + [ A{cnt4} ]'*[ A{cnt4} ];
      end;
    end;

    % keyboard;
    for cnt4=1:Nslc;                    % construct the objective matrix
      rptupd(prnt,sprintf('slc %d ',cnt4));
      At = [A{cnt4} ]';
      if (usenorm), At = At*diag(p_src{cnt4}); end;
      b2 = [];
      for l=1:ncoils,
        rptupd(prnt,'.');
        % A is the same for each coil ...
        b = vec( squeeze(Sk( plist{cnt}, kx_indx, cnt4, l )).' ); % for dual-kernel, this defaults
                                                                % to the target data
        if (usenorm)
          b = b.*p_src{cnt4};                  % normalize by the power in each linsys
        end;

        b2(:,l) = b;
        % n(:,l) = A \ b;
        % n(:,l) = Apinv * b;
        if strcmp(solver,'cgsr')
          n(:,l) = cgsolv( AA + chi*eye(size(AA)), At*b, zeros(size(At,1),1), size(At,1) );
        else
          [ n(:,l), ~ ] = bicg( AA + chi*eye(size(AA)), At*b );
        end;

        if (verbose>2), pltcmplx( At'*n(:,l), b ); keyboard; end;
      end;
      b1 = b2;
      N{cnt4}(length(N{cnt4})).n = n;
      N{cnt4}(length(N{cnt4})).chi = chi;
      N{cnt4}(length(N{cnt4})).eta = eta;

      rptupd(prnt,'\n');
    end;  
  end

end;
rptupd(prnt,'recon missing lines\n');
%% now go recon the missing lines of k-space:

if (verbose); keyboard; end;

kypts = 1:length(mss);
kypts( 1:floor(length(mss)/2) ) = kypts( (floor(length(mss)/2)):-1:1 );

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
    v = extrapolate_pattern(pattern{cnt4},nky);
%    f = identify_feedback_indices(pattern{cnt4});
%    nfc = length(f);

%    if ( (kycnt < length(mss)/2) && sum(find( pattern{cnt4} == 'm')) ), continue, end;
%    if ( (kycnt > length(mss)/2) && sum(find( pattern{cnt4} == 'p')) ), continue, end;

    indx  = cnt2 + v;
    % indx2 = cnt2 + f;
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

  if radius>1.0; radius=1.0; end;
  if radius<0.1; radius=0.1; end;
  rsz = floor(radius*size(Sk,2));

  A = zeros( rsz, length(indx3)*length(v)*size(s,3) );
  B = []; % zeros( rsz, length(indx3)*length(f)*size(s,3) );
  
  for cnt3=1:length(indx3),
    y0 = mod( indx3(cnt3) + [(size(Sk,2)-rsz)/2+(1:rsz)] , size(Sk,2) );
    y0( y0 == 0 ) = size(Sk,2);
    
    % A( (1:size(s,2)), cnt3:length(indx3):size(A,2) ) = ...
    A( (1:rsz), (cnt3-1)*nky*ncoils + (1:nky*ncoils) ) = ...
        reshape( permute( Sk(indx,y0,:), [ 2 1 3 ]), ...
                 rsz, ncoils*length(indx) );

    % B( (1:size(s,2)), (cnt3-1)*length(indx2)*ncoils + (1:length(indx2)*ncoils) ) = ...
    %     reshape( permute( Sk(indx2,y0,:), [ 2 1 3 ]), ...
    %              size(Sk,2), size(Sk,3)*length(indx2) );
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
  n = N(cnt3).n;
  
  for l=1:ncoils,
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
function v = extrapolate_pattern(p,YKS)
 
r = find( (p == 'x') | (p == '<') | (p == '>') );
v = find( (p == '+') | (p == '-') | (p == '<') | (p == '>') ) - r;

if length(v)>YKS,
  v = v( find( v~= 0 ) );
end;
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

