function [Sk,Sk2] = dpg_recon( Fin, N, dk, nky );

Sk  = zeros(size(Fin{1}));
if nargout > 1,
  Sk2 = zeros(size(Fin{1}));
end;
sz = size(Sk);


if nargin < 4,
  nky = 2;
end;

if ( ~isfield( N, 'n') & isfield( N, 'coef' ) )
  for cnt=1:length(N),
    N(cnt).n = N(cnt).coef;
  end;
end;

if isstr(nky),
  kernel = nky;
  tmp = regexp(kernel,'x');
  nky = str2num(kernel(1:(tmp-1)));
end;

nkx = size(N(1).n,1) / size(Fin{1},3) / nky;
ncoils = size(Sk,3);

% fprintf(' ky:%d kx:%d nc:%d\n',nky,nkx,ncoils);
indx3 = ceil((0-nkx)/2) -1 + (1:nkx); % [-2:2];

verbose=1;
prnt = 0;

% for kycnt=2:2*dk:(size(Sk,1)-2*dk),
for kycnt=2:2*dk:(size(Sk,1)),

  % fprintf('---%d ---\n',kycnt);
  for patcnt = 1:length(N); % [ 2:length(N)/2  (length(N)/2+(2:length(N)/2)) ];
    
    % fprintf(' %d',patcnt)
    v = extrapolate_pattern(N(patcnt).pattern,nky);

    if (patcnt <= length(N)/2)
      indx  = kycnt + (patcnt-1) + v;
    else
      indx  = kycnt + (patcnt-2) + v;
    end;        

    % fprintf(': '); fprintf(' %d', indx); fprintf('\n');
    if sum(indx>sz(1));
      indx = mod( indx, sz(1) );
      indx( indx==0 ) = sz(1);
    end;
    % fprintf('# '); fprintf(' %d', indx); fprintf('\n');

    %% recon k-space line based on that pattern:
    A = zeros( size(Sk,2), length(indx3)*length(v)*size(Sk,3) );
  
    for cnt3=1:length(indx3),
      y0 = mod( indx3(cnt3) + [1:size(Sk,2)] , size(Sk,2) );
      y0( y0 == 0 ) = size(Sk,2);

      tmpSk = zeros([ length(indx) length(y0) size(Sk,3) ]);
  
      if (patcnt <= length(N)/2)
        tmpSk(1:2:end,:,:) = Fin{1}( indx(1:2:end), y0, :);
        tmpSk(2:2:end,:,:) = Fin{2}( indx(2:2:end), y0, :);
      else
        tmpSk(1:2:end,:,:) = Fin{2}( indx(1:2:end), y0, :);
        tmpSk(2:2:end,:,:) = Fin{1}( indx(2:2:end), y0, :);
      end;

      A( (1:size(Sk,2)), (cnt3-1)*nky*ncoils + (1:nky*ncoils) ) = ...
          reshape( permute( tmpSk, [ 2 1 3 ]), ...
                   size(Sk,2), size(Sk,3)*length(indx) );

    end;

    if (verbose), % ( sum(Sk(Y(mss(cnt2))+1,:,l)) ~= 0 ); 
      tmp = kycnt+patcnt; 
      if (patcnt > length(N)/2); tmp = tmp-1; end;
      rptupd(prnt,sprintf(' %3d',[ tmp; indx(:) ]));
      rptupd(prnt,sprintf(' (pattern: %s)', N(patcnt).pattern) );
      rptupd(prnt,'\n');
      if (verbose>1),
        figure(1); imagesc( log10(abs(Sk(:,:,1)) + 1e-1) ); pause(0.5);
        keyboard;
      end;
    end;
    % linerpt(cnt2) = 'x';

    n = N(patcnt).n;
    
    for l=1:size(Sk,3),
      % for each coil ...

      % if cnt2==109, keyboard; end;
      % rptupd(prnt,(' line %3d, pattern %3d\n', cnt2, plist(cnt3) );

      if (nargout == 1),
        tmp = kycnt+(patcnt);
        if (patcnt > length(N)/2), tmp = tmp-1; end;
        Sk( tmp, :, l ) = [A] * n(:,l);
      else
        if (patcnt > length(N)/2 ), 
          Sk2( kycnt+(patcnt-1), :, l ) = [A] * n(:,l);
        else
          Sk( kycnt+(patcnt), :, l ) = [A] * n(:,l);
        end;
      end;
    end;

    % if (verbose>1), keyboard; end;
  end;

end;

Sk = Sk(1:sz(1),:,:);

%% need to fix this above, but for now...
Sk = Sk([2:sz(1) 1],:,:);               % rotate the output back 1 kspace line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rptupd( prnt, str )

if (prnt)
  fprintf( str );
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = extrapolate_pattern(p,YKS)
 
r = find( (p == 'x') | (p == '<' ) | (p == '>') );
v = find( (p == '*') | (p == '<' ) | (p == '>') ) - r;

if length(v)>YKS,
  v = v( find( v~= 0 ) );
end;
return;
