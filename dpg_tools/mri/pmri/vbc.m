% function b0 = vbc( Fk, mode )
%
% combines multi-coil data while preserving phase and giving more
% homogenous illumination
% 
% Fk - input k-space.  Can be either MxNxL or MxNxZxLxT...
% 
% an implementation of - Martin Buehrer, Peter Boesiger, and Sebastian Kozerke. Virtual body coil calibration for phased-array imaging. In Proceedings 17th ISMRM Scientific Meeting, page 759, Honolulu, HI, 2009.
% 

function [b0,sv] = vbc( Fk, mode )

if (nargin<2),
  mode = 'nofft';
end;

sz = size(Fk);
if length(sz) > 3,
  
  [m n z c t] = size(Fk);

  b0 = zeros([ m n z 1 t]);
  sv = zeros([ c t ]);
  for cntT = 1:t,
    [b0(:,:,:,:,cntT),sv(:,cntT)] = vbc_set( squeeze( Fk(:,:,:,:,cntT) ), mode );
  end;
  
  sz = size(Fk);
  if length(sz) > 4,
    b0 = reshape( b0, [ m n z 1 sz(5:end) ]);
    sv = reshape( sv, [ c sz(5:end) ]);
  end;
else
  [b0,sv] = vbc_set( Fk, mode );
end;

return;

function [b0,sv] = vbc_set( Fk, mode )


if length(size(Fk)) == 3,
  [m n c ] = size(Fk);
  z = 1;
else,
  [m n z c ] = size(Fk);
end;

if ~strcmp(mode,'nofft')
  Fk = fftshift(ifft2(fftshift(Fk)));
end;

Fk2 = reshape( Fk, [m*n*z c] );

d = Fk2'*Fk2;
% Fk3 = orth( Fk2 * diag(1./sqrt(diag(d))) );
% Fk3 = lsv( Fk2 * diag(1./sqrt(diag(d))) );
  
[u s v] = svd( Fk2 * diag(1./sqrt(diag(d))), 0 );
s = diag(s);

if ( (s(1)-s(2))/s(1) < 0.05 ),
  u1 = phzshift( reshape(u(:,1),[m n z]), ...
                 reshape(u(:,2),[m n z]), {'nofft'} ) / 2;
  v1 = (v(:,1)+v(:,2))/2;
       % phzshift( reshape(v(:,1),[m n z]), ...
       %           reshape(v(:,2),[m n z]), {'nofft'} ) / 2;

  s1 = (s(1)+s(2))/2;
  %  keyboard;
else
  u1 = u(:,1);
  v1 = v(:,1);
  s1 = s(1);
end;

mag = norm(vec( sqrt(sum(abs(Fk).^2,3)) ));

b0 = mag * reshape( u1, [ m n z ] );


% sv = % [Fk2 * diag(1./sqrt(diag(d))) ]' * y;
sv = mag * diag(1./sqrt(diag(d)))*v1*(1/s1);

return;
