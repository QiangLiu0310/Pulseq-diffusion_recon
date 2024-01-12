function varargout = grog( varargin )
% grog: GRAPPA Operator Gridding (for PROPELLER)
%
% based on the description given in
%
%  Non-Cartesian data reconstruction using GRAPPA operator gridding (GROG)
%  by N. Seiberlich, F. A. Breuer, M. Blaimer, K. Barkauskas, P. M. Jakob and M. A. Griswold. 
%  Magn Reson Med.  58(6):1257-1265,  Dec 2007
%
%
%
% The foundation of the method is that grappa can be used to formulate new
% points.  First, for an typical 1-delta-k shift with a 1x1 kernel, the
% sythesis equation is
%       s( kx, ky + 1 ) = G s( kx, ky ), where 
%  s : vector of data for multiple coils (size: L -by- 1)
%  G : grappa operator (size: L -by- L )
%
% Second, G can be decomposed via svd or eigen-decomposition, and then scaled.
%   G = E V E'.   G_d = E V^{1/d} E'
%
% this allows arbitrary points to be generated via
%       s( kx, ky + d ) = G_d s( kx, ky )
%
% Further, grappa operators are seperable, so one can apply combinations of
% G along x and y for regridding.
%


k = varargin{1};                          % k-space data samples, N points -by- L coils
r = varargin{2};                          % [x,y,z] locations for each k-space data point, N points -by- ndim
rnew = [];
G = [];
F = [];
verbose = 0;
kernel = '1x1'; 
refang = 0;

chi = 1e-6;

if length(varargin) < 2,
  error('need to declare both input kspace data and data point locations\n');
end;
if length(varargin) > 2,
  G = varargin{3};
end;

[ npts, ncoils ] = size(k);
[ npts1, ndim ] = size(r);
if npts~=npts1,
  error('number of k-space data points (%d) and number of data locations (%d) need to match',n,n1);
end;

vararg = varargin(4:end);
for cnt=1:2:length(vararg),
  if strcmp( vararg{cnt}, 'v' ),
    verbose = cell2mat(vararg(cnt+1));
  elseif strcmp( vararg{cnt}, 'N' ),
    G = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'kernel' ),
    kernel = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'refang' ),
    refang = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'acs' ),
    acs = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'rnew' ),
    rnew = vararg{cnt+1};
  elseif strcmp( vararg{cnt}, 'chi' ),
    chi = vararg{cnt+1};
  end;
end;


if isempty(G),
  % compute the grappa gridding coefficients

%   if exist('cgsolv') ~= 2,
%     %% add a path to the conjugate gradient solver, if needed.
%     ans = which('grog');
%     idx = findstr('mri/parallel',ans)-1;
%     addpath([ ans(1:idx) 'mathlib/conjgrad' ]);
%   end;
  
  d = double(int32(r)) - r;
  if (sum(d(:))~=0)
    error('for calibration stage, all input data locations must be integers\n');
  end;
  for l=1:ncoils,
    kcal(:,:,l) = full(sparse( r(:,1),r(:,2),k(:,l) )); % ,max(r(:,2)),max(r(:,2)) ));
  end;

  nky = str2double(kernel(1));          % nky : number of points in kernel along y direction
  nkx = str2double(kernel(3:end));      % should be odd
  indx3 = ceil((0-nkx)/2) -1 + (1:nkx); % [-2:2];

  if (1),
    tmp = repmat('x',1,nky);

    % pattern{1} = [ tmp 'o*' ];
    pattern{1} = [ tmp '*'  ];
    pattern{2} = [  '*' tmp ];
    % pattern{4} = [ '*o' tmp ];
  end;
  
  kcal2 = sum( kcal~=0, 3 );
  kcal2( kcal2 ~= ncoils ) = 0;
  kcal2( kcal2 == ncoils ) = 1;
  
  Y{1} = find( sum(kcal2,1) );          % x sample locations
  Y{2} = find( sum(kcal2,2) );          % y sample locations
  
  acsx{1} = compute_acs(Y{1});
  acsx{2} = compute_acs(Y{2});

  if 1, % for dircnt=2:-1:1,

%     if dircnt==1,
%       Sk = permute(kcal(acsx{2},:,:),[2 1 3]);
%       % Sk = permute(kcal,[2 1 3]);
%     else,
%      Sk = kcal(:,acsx{1},:);
      % Sk = kcal;
%    end
    s = kcal; % (Y{2},acsx{1},:);
    
    acs = acsx{2};
    
    for patcnt = 1:length(pattern),
      v = extrapolate_pattern( pattern{patcnt} );
      plist = [];
      
      for cnt2 = 1:length(acs);
        %% check if the GRAPPA pattern...
        indx = acs(cnt2) + v;
        
        %% can be accomodated by the sampling pattern 
        tmp = find( Y{2}(:)*ones(1,length(indx)) == ones(length(Y{2}),1)*indx );

        if ( length(tmp) == length(indx) ),
          plist =  [ plist acs(cnt2) ];
        end;
      end;
      
      A = zeros( length(plist)*size(s,2), length(indx3)*nky*ncoils );
      n = zeros( (nky)*ncoils*length(indx3), ncoils*nky*length(indx3) );

      for cnt2 = 1:length(plist),
        %% and ACS line ...

        %% if the sampling pattern accomodates the GRAPPA pattern...
        indx  = plist(cnt2) + v;         % feedforward indices

        %% determine the k-space filling coefficients 
        for cnt3=1:length(indx3),
          y0 = mod( indx3(cnt3) + [1:size(s,2)] , size(s,2) );
          y0( y0 == 0 ) = size(s,2);
          
          % A( size(s,2)*(cnt2-1) + (1:size(s,2)), cnt3:length(indx3):size(A,2) ) = ...
          A( size(s,2)*(cnt2-1) + (1:size(s,2)), (cnt3-1)*nky*ncoils+(1:nky*ncoils) ) = ...
            reshape( permute( s(indx,y0,:), [ 2 1 3 ]), size(s,2), size(s,3)*length(indx) );
        end;
        %% for an N-by-1 kernel, these next lines are enough...
        % A2 = [ A2; ...
        %         reshape( permute(squeeze( Sk(indx,:,:) ),[ 2 3 1 ]), ...
        %                  size(s,2), size(s,3)*length(indx) )
        %       ];
        % 
        % b2 = [b2; squeeze( Sk( plist{cnt}(cnt2), :, l ) ).' ];
      end;

      %% permute A to match the operator mapping to I for d=0
      % idx = vec( reshape( 1:size(A,2), nkx, ncoils ).' );
      % A = A(:,idx);
      
      % Apinv = pinv(A);
      % AA = A'*A;
      AA = [A]'*[A];
      At = [A]';

      for dircnt=2:-1:1,
        if ((dircnt==1) && (dx~=0) ),
          plist2 = mod( plist + v, size(s,1) );
          plist = plist2; clear plist2;
          plist( plist == 0 ) = size(s,1);
          dx = -v;
        else
          dx = 0;
        end;
        for cnt3=1:nkx, % nkx = length(indx3),
          for l=1:ncoils,
            % rptupd(prnt,'.');
            y0 = mod( (1:size(s,2))+indx3(cnt3) + dx, size(s,2) );
            y0( y0 == 0 ) = size(s,2);
            % A is the same for each coil and target point ...
            b = vec( squeeze(s( plist, y0, l )).' );
        
            % n(:,l) = A \ b;
            % n(:,l) = Apinv * b;
            % n( :, (cnt3-1)*ncoils + l ) = cgsolv( AA, At*b, zeros(size(At,1),1), size(At,1) );
            n( :, (cnt3-1)*ncoils + l ) = bicg( AA + chi*eye(size(AA)), At*b );
            
            if (verbose>2), pltcmplx( At'*n(:,l), b ); keyboard; end;
          end;
        end;
        
        G(length(G)+1).pattern = pattern{patcnt};
        if (dircnt == 2),
          G(end).dir = 'normal';
        else, 
          G(end).dir = 'parallel';
        end;
        G(end).n   = n;         % Jan 2011: the 'sign' of the coefficients depends on
                                % whether k-space is 'chopped' or not.  raw
                                % k-space makes an image with ffts(ifft2(ffts(
                                % .. ))).  The FIL demo data had its phase
                                % altered to remove the need for the last ffts.
                                % For some reason, this affects the GROG
                                % coefficients here, necessitating a minus sign
                                % on the coefs. OTW, the eigen-value scaling
                                % doesn't work.
        
        G(end).refang = refang;              % (x,y) offset for reference data angle (should be in radian, [0 .. pi/2])
        G(end).kernel = kernel;
      end;
    end;
  end;
  % end of grog calibration section  

  varargout{1} = G;
  
else,
  % perform the gridding



  r_orig = r;
  if isempty(rnew)
    rnew = floor(r);
  else
  end;
  % rnew = floor(r+0.5)-0.5;   %round? fix?         % points on Cartesian grid

  knew = zeros([ size(rnew,1) size(k,2) ]);
  
  if(0)
  %%% patch any holes that remain in the target grid
  mn = min( rnew ) - 1;
  rnew2 = rnew - ones(size(rnew,1),1)*mn;
  Z     = zeros( round(max(rnew2(:,1))), round(max(rnew2(:,2))) );

  for cnt=1:size(rnew2,1),
    Z(rnew2(cnt,1),rnew2(cnt,2)) = 1;
  end;
  
  intpt=zeros(size(Z)); 
  for scnt=2:(length(Z(:))-1); 
    if (Z(scnt)==0)&(Z(scnt-1)~=0)&(Z(scnt+1)~=0); 
      intpt(scnt)=1; 
    end; 
  end;
  [xx,yy] = find(intpt);

  if length(xx)>0,
    for cnt=1:length(xx);
      tmp = [xx(cnt) yy(cnt)] + mn;
      for cnt2=1:size(r,1); 
        d(cnt2) = sqrt(sum( abs([tmp-r(cnt2,:)]).^2 )); 
      end;
      [ttt,iii]=min(d);
      snew(cnt,:) = tmp;
      s(cnt,:) = r(iii,:);
    end;
  
    rnew = [ rnew; snew ];
    r    = [ r; s ];
  end;
  % fprintf('length r: %d %d',length(r_orig),length(r));
  end;

  % dist = rnew - r;                      % distance metric to move points
  % sz = round( max(r) - min(r) );

  for cnt=1:length(G),
    v(cnt,:) = extrapolate_pattern(G(cnt).pattern);
    hrz(cnt) = strcmp( G(cnt).dir, 'parallel' );
    vrt(cnt) = strcmp( G(cnt).dir, 'normal' );
  
    [ V D ] = eig( G(cnt).n );
    G(cnt).e = D;
    G(cnt).ev = V;
  end;
  
  for cnt=1:size(rnew,1), % npts, % npts/2+32, % 

    % if( (rnew(cnt,1)==0) & (rnew(cnt,2)==0) ), keyboard; end;

    %% project distance along the reference angle directions
    
    %%     find the closest source point:
    [tmp,srcindx]=min( sqrt(sum((ones(length(r),1)*rnew(cnt,:) - r).^2,2)));
    dist = rnew(cnt,:) - r(srcindx,:);
    
    dp = dist * [ sin(refang) cos(refang) ]'; % distance parallel
    dn = dist * [ cos(refang) -sin(refang) ]'; % distance normal
    
    %% compute the operator matrix, combining the different dimmensions
    %% into a single operator
    O = 1;                              % initialize the grappa operator
    for dircnt=2:-1:1,
      coef = [];
      if dircnt==1,
        if ( dn > 0 )  % vertical disp
          coef = find( ( sum(v == -1,2)~=0 ) & (vrt(:) == 1) );
        elseif ( dn < 0 )  % vertical disp
          coef = find( ( sum(v ==  1,2)~=0 ) & (vrt(:) == 1) );
        end;
        ddd = dn;
      else,
        if ( dp > 0 )  % horizontal disp
          coef = find( ( sum(v == -1,2)~=0 ) & (hrz(:) == 1) );
        elseif ( dp < 0 )  % horizontal disp
          coef = find( ( sum(v ==  1,2)~=0 ) & (hrz(:) == 1) );
        end;
        ddd = dp; 
      end;
%      fprintf(' %d %f\n',coef,ddd);
      
      if ~isempty(coef)
        O = O * ([G(coef).ev * diag(diag(G(coef).e).^abs( ddd )) * inv( G(coef).ev )]);
        kernel = G(coef).kernel;
        refang = G(coef).refang;
      end;
    end;

    nky = str2double(kernel(1));           
    % indx  = v;
    nkx = str2double( kernel(3:end) ); % should be odd
    indx3 = ceil((0-nkx)/2) -1 + (1:nkx); % [-2:2];
     
%     %% recon k-space line based on that pattern:
%     A = zeros( size(s,2), length(indx3)*length(v)*size(s,3) );
%     B = zeros( size(s,2), length(indx3)*length(f)*size(s,3) );
%   
%     for cnt3=1:length(indx3),
%       y0 = mod( indx3(cnt3) + [1:size(Sk,2)] , size(Sk,2) );
%       y0( y0 == 0 ) = size(Sk,2);
%     
%       A( (1:size(s,2)), cnt3:length(indx3):size(A,2) ) = ...
%           reshape( permute( Sk(indx,y0,:), [ 2 1 3 ]), ...
%                    size(Sk,2), size(Sk,3)*length(indx) );
%     end;

    % compute O_y
    % keyboard;
%    idx = vec( reshape( 1:ncoils*nkx, ncoils, nkx ).' );


    %%% 1st attempt:
    % 1. find all points that are |1| or |2| delta-k away from r(cnt,:)
    % 2. of those points, find which are colinear (phase is equal-n-opposite) with reference/calibration angle

    % dlt = sqrt(sum( abs( [ r - ones( size(r,1), 1 )*r(cnt,:) ] ).^2 , 2 ));
    % idx = find( dlt == abs(indx3(cnt2)) );
    
    %%% 2nd (and 0th) attempt:
    % 1. calculate the reference points needed, from the indx3 and refangle variables
    % 2. find those points w/in the r list
    
    tol = 1e-8;
    if length(O)~=1,
      
      % compute coordinates of needed reference points
      y0 = indx3' * [ sin(refang) cos(refang) ] + ones(nkx,1)*r(srcindx,:);
      y1 = y0;
    
      idx2 = y0(:,2) < min(r(:,2)); 
      y0(idx2,2) = y0(idx2,2) + ( max(r(:,2))-min(r(:,2)) + 1 );

      idx2 = (y0(:,2) - max(r(:,2))) > tol;
      y0(idx2,2) = y0(idx2,2) - ( max(r(:,2))-min(r(:,2)) + 1 );

      idx1 = y0(:,1) < min(r(:,1));
      y0(idx1,1) = y0(idx1,1) + ( max(r(:,1))-min(r(:,1)) + 1 );

      idx1 = (y0(:,1) - max(r(:,1)) ) > tol;
      y0(idx1,1) = y0(idx1,1) - ( max(r(:,1))-min(r(:,1)) + 1 );

      srcpts = [];
      for cnt2 = 1:length(indx3),
        tmp = find(  sum( abs( r_orig - ones(size(r_orig,1),1)*y0(cnt2,:) ) < tol, 2) == 2 ); % find that point in the reference data
        if ~isempty(tmp)
          srcpts = [ srcpts tmp ];        % if found, add this point to the list.
        end;
      end;
      % apply the operator to the data
      % disp(srcpts);
      if length(srcpts) == nkx,
        knew(cnt,:) = vec( k( srcpts,:).' ).' * O( :, ((find(indx3==0)-1)*ncoils)+(1:ncoils) );
      else,
        knew(cnt,:) = 0;
        %      keyboard;
      end;
      
    else,
      knew(cnt,:) = k(cnt,:);
    end;
    
  end;

  % Now, populate the Cartesian matrix
  mn = min( rnew ) - 1;

  rnew2 = rnew - ones(size(rnew,1),1)*mn;
  
  F    = zeros( round(max(rnew2(:,1))), round(max(rnew2(:,2))), ncoils );
  scl  = zeros([ size(F,1) size(F,2) ]); 

  for cnt=1:size(rnew2,1),
    
    F( rnew2(cnt,1), rnew2(cnt,2), : ) = ...
        reshape( F( rnew2(cnt,1), rnew2(cnt,2), :), [1 ncoils] ) + knew(cnt,:);
    if sum( knew(cnt,:)~=0 ) == ncoils,
      scl( rnew2(cnt,1), rnew2(cnt,2) ) = scl( rnew2(cnt,1), rnew2(cnt,2) ) + 1;
    end;
  end;

  scl2 = scl;
  scl2(scl==0) = 1;
  for lcnt=1:ncoils,
    F(:,:,lcnt) = F(:,:,lcnt) ./ scl2;
  end;

  if nargout>0, varargout{1} = knew; end;
  if nargout>1, varargout{2} = rnew; end;
  if nargout>2, varargout{3} = F;    end;
  
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = extrapolate_pattern(p)
 
r = find( p == 'x' );
v = find( p == '*' ) - r;
return;

function acs = compute_acs( Y )

d = diff( Y(:) );

acs = Y([ min(find(d==1)); find( d == 1)+1 ]);

