function [F,I,N] = dpg_cal( varargin );

% [F,I,N ] = dpg_cal( Wk, Sk, Y, [opt args]);
% 
% opt args:
%  'v', {0,1}              - verbosity off/on
%  'kernel', string        - grappa kernel size, x by y. e.g. '2x5', '4x1'
%  'N', (xy)-by-num_patterns-L matrix   - grappa reconstruction parameters
%  'acs', array of indices - index locations to use for ACS lines
% 
% This implementation is derived from 
% "Generalized autocalibrating partially parallel acquisitions"
%   by M. A. Griswold, P. M. Jakob, et. al..
%  Mag. Reson. Med, 47(6):1202-1210, Jun 2002
%

%
% Copyright 2004-2016  Scott Hoge  (shoge at ieee dot org)
%
% All rights reserved. Licensed according to the GNU General Public Licence
% (see http://www.gnu.org)
%

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

Wk  = varargin{1};
kin = varargin{2};
Y   = varargin{3}; % the size of PE dimension

acs = []; N = []; verbose=0; kernel = '2x5'; dks = []; prnt = 0;
eta = 1.0;
chi = 0.0001;

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
  elseif (strcmp( vararg{cnt}, 'report' ) || strcmp( vararg{cnt}, 'r' )),
    prnt = vararg{cnt+1};
  end;
end;

if length(size(Wk)) == 3, 
  [Wm,Wn,ncoils] = size(Wk);
else,
  Wm = Wk(1); Wn = Wk(2); ncoils = Wk(3);
end;


%% zero-pad the k-space data
Sk = zeros( Wm, Wn, ncoils );
if iscell(kin)
  dualks = 1; % two k-space kernels, lq
  s = kin{3};
else
  dualks = 0;
  s = kin; clear kin;
end;
Sk(Y,:,:) = s;

rptupd(prnt,['* ' mfilename ': kernel size = ' kernel '\n']);

d = diff(Y(:));
if ( isempty(acs) ),
  prfl = max( sum(abs( Sk(:,:,:)),3).' ); % max x for each lines in averages slices? lq

  if length( d == 1) == 1, 
    acs = find( prfl == max(prfl(Y(find(d==1)))) ); 
  else
    acs = Y([ min(find(d==1)); find( d == 1)+1 ]); % this line output: Y[1;2], lq
  end;
  rptupd(prnt,['* ' mfilename ': acs =']);
  rptupd(prnt,sprintf(' %d', acs) );
  rptupd(prnt,'\n');
end;

if ( isempty(dks) ), % if (1), nevermind, lq
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
  for cnt=1:length(N)
    pattern{cnt} = N(cnt).pattern;
  end;
end;

if isempty(pattern)
%% foreach dk spacing
for cnt=1:length(dks),
  
  if ( dks(cnt) == 1 )
    % for no acceleration, just specify the pattern
    tmppat = repmat('*',[1 nky]);
    for patcnt = 1:nky,
      pattern{length(pattern)+1} = tmppat;
      pattern{length(pattern)}(patcnt) = '<';
    end;

    for patcnt = 1:nky,
      pattern{length(pattern)+1} = tmppat;
      pattern{length(pattern)}(patcnt) = '>';
    end;

    continue;
  end;

  tmp = [ repmat([ '*' repmat('o',1,dks(cnt)-1) ],1,nky-1) '*' ];
  z2 = strfind(tmp,'o');
  if isempty(z2),                       
    % pad with zeros, if the kernel is 1-by-N
    tmp = [ 'o' tmp 'o' ];   
    z2 = strfind(tmp,'o');
  end;
  [tmp2,indx] = sort( abs(z2 - z2( floor(median(1:length(z2))) )) );

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
      if ( cnt2 <= npat/2 )
        pattern{length(pattern)}( cnt3 ) = '<';
      else
        pattern{length(pattern)}( cnt3 ) = '>';
      end;
    end;
  
  end;

end;
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

if isempty(N),
  % N = zeros( length(indx3)*nky*size(s,3),length(pattern),size(s,3));

  rptupd(prnt,'calc recon param\n');
  for cnt = 1:length(pattern), 
    v = extrapolate_pattern(pattern{cnt},nky);

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
    A = zeros( length(plist{cnt})*size(s,2), length(indx3)*nky*size(s,3) );
    
    rptupd(prnt,sprintf('%3d  %s: ', cnt, pattern{cnt}));

    % v = extrapolate_pattern(pattern{cnt});
    N(length(N)+1).pattern = pattern{cnt};

    n = zeros( nky*size(s,3)*length(indx3),size(s,3));

    for cnt2 = 1:length(plist{cnt}),
      %% and ACS line ...

      %% if the sampling pattern accomodates the GRAPPA pattern...
      indx  = plist{cnt}(cnt2) + v;         % feedforward indices
      % indx2 = plist{cnt}(cnt2) + f;         % feedback indices

      %% determine the k-space filling coefficients 
      for cnt3=1:length(indx3),
        y0 = mod( indx3(cnt3) + [1:size(Sk,2)] , size(Sk,2) );
        y0( y0 == 0 ) = size(Sk,2);
        
        % generate a hybrid source matrix
        if (dualks)
          tmpSksub = zeros( size( Sk(indx,y0,:) ) );
          if ( sum(pattern{cnt}=='<') ), 
            tmpSksub( 1:2:length(indx), :, : ) = kin{1}( indx(1:2:end),y0,: );
            tmpSksub( 2:2:length(indx), :, : ) = kin{2}( indx(2:2:end),y0,: );
%            tmpSksub( (v < 0), :, : )  = kin{2}( indx( (v < 0) ),y0,: );
%            tmpSksub( (v >= 0), :, : ) = kin{3}( indx( (v >= 0) ),y0,: );
          else
            tmpSksub( 1:2:length(indx), :, : ) = kin{2}( indx(1:2:end),y0,: );
            tmpSksub( 2:2:length(indx), :, : ) = kin{1}( indx(2:2:end),y0,: );
%            tmpSksub( (v >= 0), :, : ) = kin{2}( indx( (v >= 0) ),y0,: );
%            tmpSksub( (v < 0), :, : )  = kin{3}( indx( (v < 0) ),y0,: );
          end
        else
          tmpSksub = Sk(indx,y0,:);
        end;
        % A( size(s,2)*(cnt2-1) + (1:size(s,2)), cnt3:length(indx3):size(A,2) ) = ...
        A( size(s,2)*(cnt2-1) + (1:size(s,2)), (cnt3-1)*nky*ncoils + (1:nky*ncoils) ) = ...
            reshape( permute( tmpSksub, [ 2 1 3 ]), size(Sk,2), size(Sk,3)*length(indx) );

      end;
      %% for an N-by-1 kernel, these next lines are enough...
      % A2 = [ A2; ...
      %         reshape( permute(squeeze( Sk(indx,:,:) ),[ 2 3 1 ]), ...
      %                  size(s,2), size(s,3)*length(indx) )
      %       ];
      % 
      % b2 = [b2; squeeze( Sk( plist{cnt}(cnt2), :, l ) ).' ];
    end;

    % normalize, following to Algo 4.2 of M. Schneider, "GPGPU for Accelerated GRAPPA Autocalibration in
    %    Magnetic Resonance Imaging," Masters Thesis, U of Erlangen, Apr 2008
    % A_bak = A;
    p_src = sum(abs(A),2).^(-eta/2);    % compute the power in each line of the linsys
    A = tmult(A,diag(p_src),1);         % normalize by the power in each linsys line

    % Apinv = pinv(A);
    % AA = A'*A;
    AA = [A ]'*[A ];
    At = [A ]';
    clear b2;
    for l=1:size(Sk,3),
      rptupd(prnt,'.');
      % A is the same for each coil ...
      b = vec( squeeze(Sk( plist{cnt}, :, l )).' ); % for dual-kernel, this defaults
                                                    % to the target data

      b = b.*p_src(:);                  % normalize by the power in each linsys line

      b2(:,l) = b;
      n(:,l) = A \ b;

      if (verbose>2), pltcmplx( At'*n(:,l), b ); keyboard; end;
    end;
    b1 = b2;
    N(length(N)).n = n;
    
    rptupd(prnt,'\n');
  end;  
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

    if ( (kycnt < length(mss)/2) && sum(find( pattern{cnt4} == '-')) ), continue, end;
    if ( (kycnt > length(mss)/2) && sum(find( pattern{cnt4} == '+')) ), continue, end;

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
  rsz = size(Sk,2);
  A = zeros( rsz, length(indx3)*length(v)*size(s,3) );
  B = []; % zeros( rsz, length(indx3)*length(f)*size(s,3) );
  
  for cnt3=1:length(indx3),
    y0 = mod( indx3(cnt3) + [(size(Sk,2)-rsz)/2+(1:rsz)] , size(Sk,2) );
    y0( y0 == 0 ) = size(Sk,2);
    
    A( (1:rsz), (cnt3-1)*nky*ncoils + (1:nky*ncoils) ) = ...
        reshape( permute( Sk(indx,y0,:), [ 2 1 3 ]), ...
                 rsz, size(Sk,3)*length(indx) );
  end;

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
function v = extrapolate_pattern(p,YKS)
 
r = find( (p == 'x') | (p == '<') | (p == '>') );
v = find( (p == '*') | (p == '<') | (p == '>') ) - r;

if length(v)>YKS,
  v = v( find( v~= 0 ) );
end;
return;


