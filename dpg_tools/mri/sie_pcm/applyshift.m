function b = applyshift(a,x,y)

[m n c]=size(a);

b = zeros(size(a));
for cnt=1:c,
  A = fftshift(fft2(fftshift( a(:,:,cnt) )));
  
  A = diag( exp( j*2*pi*x/m*[(1:m) - m/2]) ) * A;
  A = A * diag( exp( -j*2*pi*y/n*[(1:n) - n/2]) );
  
  b(:,:,cnt) = ifftshift( ifft2( ifftshift( A ) ) );
end;

