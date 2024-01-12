function [Sk,Sk2] = dpg_recon( Fin, N, dk, nky );

if iscell(Fin);
  ROp = Fin{1};
  ROn = Fin{2};
elseif isstruct(Fin);
  ROp = Fin.p;
  ROn = Fin.n;
else
  error('first input (Fin) needs to be a structure with Fin.p and Fin.n');
end;
clear Fin;

if ( ~isfield( N, 'n') & isfield( N, 'coef' ) )
  for cnt=1:length(N),
    N(cnt).n = N(cnt).coef;
  end;
end;

nseg = 1;
if iscell(ROp)
  % for multi-source segEPI images, data from each segment are stored in a
  % different cell.
  sz = size(ROp{1});
  nseg = length(ROp);
else
  sz = size(ROp);
end;
Sk  = zeros([ sz size(N(1).n,3) ]);
if nargout > 1,
  Sk2 = zeros(size(Sk));
end;



if nargin < 4,
  nky = 2;
end;

if isstr(nky),
  kernel = nky;
  tmp = regexp(kernel,'x');
  nky = str2num(kernel(1:(tmp-1)));
end;

nkx = size(N(1).n,1) / sz(3) / nky;
ncoils = size(Sk,3);

% seg = nky / dk;

% fprintf(' ky:%d kx:%d nc:%d\n',nky,nkx,ncoils);
indx3 = ceil((0-nkx)/2) -1 + (1:nkx); % [-2:2];

verbose=0;
prnt = 0;

bins = zeros([ 1 size(Sk,1) ]);
if iscell(ROp)
  for cnt=1:length(ROp)
    bins( find( sum(sum(abs(ROp{cnt}),3),2) ) ) = +cnt;
    bins( find( sum(sum(abs(ROn{cnt}),3),2) ) ) = -cnt;
  end;
else
  bins( find( sum(sum(abs(ROp),3),2) ) ) = +1;
  bins( find( sum(sum(abs(ROn),3),2) ) ) = -1;
end;
% keyboard;
% for kycnt=2:2*dk:(size(Sk,1)-2*dk),
% for kycnt=2:2*dk:(size(Sk,1)),
for kycnt=1:size(Sk,1),

  if (prnt) fprintf('---%d ---\n',kycnt); end;
  for patcnt = 1:length(N); % [ 2:length(N)/2  (length(N)/2+(2:length(N)/2)) ];
    
    if (prnt) fprintf(' %d',patcnt); end;
    v = extrapolate_pattern(N(patcnt).pattern,nky);

    ptrn = N(patcnt).pattern;

    p = find( (ptrn=='+') | (ptrn=='-') | (ptrn=='<') | (ptrn=='>') );

    indxpos = find( (ptrn(p) == '+') | (ptrn(p) == '>') );
    indxneg = find( (ptrn(p) == '-') | (ptrn(p) == '<') );
    
    indx  = kycnt + v;

    % fprintf(': '); fprintf(' %d', indx); fprintf('\n');
    if sum(indx>sz(1));
      indx = mod( indx, sz(1) );
    end;
    if sum(indx<0);
      indx = mod( indx+sz(1), sz(1) );
    end;
    indx( indx==0 ) = sz(1);
    % fprintf('# '); fprintf(' %d', indx); fprintf('\n');

    % if (kycnt == size(Sk,1)); fprintf(' %s',ptrn); end;
    % if ( (kycnt == size(Sk,1)) && strcmp('>oo-',ptrn) ); keyboard; end;
    
    if ( sum([ (bins(indx(indxpos)) == 1:nseg) ...
               ( -bins(indx(indxneg)) == 1:nseg) ]) ~= length(indx) )
      continue;
    end;
         


    %% recon k-space line based on that pattern:
    A = zeros( size(Sk,2), length(indx3)*length(v)*size(Sk,3) );
  
    for cnt3=1:length(indx3),
      y0 = mod( indx3(cnt3) + [1:size(Sk,2)] , size(Sk,2) );
      y0( y0 == 0 ) = size(Sk,2);

      tmpSk = zeros([ length(indx) length(y0) size(Sk,3) ]);
  
      if iscell(ROp)
        for cnt4=1:nseg
          tmpSk(indxpos(cnt4),:,:) = ROp{cnt4}( indx(indxpos(cnt4)), y0, :);
          tmpSk(indxneg(cnt4),:,:) = ROn{cnt4}( indx(indxneg(cnt4)), y0, :);
        end;
      else,
        tmpSk(indxpos,:,:) = ROp( indx(indxpos), y0, :);
        tmpSk(indxneg,:,:) = ROn( indx(indxneg), y0, :);
      end;
      A( (1:size(Sk,2)), (cnt3-1)*nky*ncoils + (1:nky*ncoils) ) = ...
          reshape( permute( tmpSk, [ 2 1 3 ]), ...
                   size(Sk,2), size(Sk,3)*length(indx) );

    end;

    if (verbose), % ( sum(Sk(Y(mss(cnt2))+1,:,l)) ~= 0 ); 
      tmp = kycnt+patcnt; 
      % if (patcnt > length(N)/2); tmp = tmp-1; end;
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
        tmp = kycnt; % +(patcnt);
        % if (patcnt > length(N)/2), tmp = tmp-1; end;
        Sk( tmp, :, l, : ) = tmult( n(:,l,:), A, 1 );
      else
        if (patcnt > length(N)/2 ), 
          Sk2( kycnt+(patcnt-1), :, l, : ) = tmult( n(:,l,:), A, 1 );
        else
          Sk( kycnt+(patcnt), :, l, : ) = tmult( n(:,l,:), A, 1 );
        end;
      end;
    end;

    % if (verbose>1), keyboard; end;
    break;
  end;

end;

Sk = Sk(1:sz(1),:,:,:);
%Sk  = Sk([2:sz(1) 1],:,:,:);
%if nargout > 1,
%Sk2 = Sk2([2:sz(1) 1],:,:,:);
%end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rptupd( prnt, str )

if (prnt)
  fprintf( str );
end;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = extrapolate_pattern(p,YKS)
 
r = find( (p == 'x') | (p == '<' ) | (p == '>') );
v = find( (p == '*') | (p == '+') | (p == '-') | (p == '<' ) | (p == '>') ) - r;

if length(v)>YKS,
  v = v( find( v~= 0 ) );
end;
return;
