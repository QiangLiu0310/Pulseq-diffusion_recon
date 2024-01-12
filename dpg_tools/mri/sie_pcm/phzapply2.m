function [a_out,b_out] = phzapply2(a,b,y);
% [a_out] = phzapply2(a,b,y);
%
% apply the EPI phase shift parameters calcuated by phzshift to 
% seperate RO+ and RO- data
%


[m n c]=size(a);

[pm pc ] = size(y);

A1 = a; % a(1:2:end,:,1:c);
A1a = ffts(ifft(ffts(A1,[ 2]),[],2),[ 2]);

A2 = b; % a(2:2:end,:,1:c);
A2a = ffts(ifft(ffts(A2,[ 2]),[],2),[ 2]);

A1b = zeros(size(A1a));
A2b = zeros(size(A2a));
for cnt=1:c,
  A1b(:,:,cnt) = tmult( A1a(:,:,cnt), diag( exp( +j*y(1,cnt)/2*((1:size(A2a,2))-size(A2a,2)/2))), 2 );
  A2b(:,:,cnt) = tmult( A2a(:,:,cnt), diag( exp( -j*y(1,cnt)/2*((1:size(A2a,2))-size(A2a,2)/2))), 2 );
end;

A1c = zeros(size(A1a));
A2c = zeros(size(A2a));
for cnt=1:c,
  A2c(:,:,cnt) = iffts(fft(iffts(A2b(:,:,cnt),[ 2]),[],2),[ 2]);
  A1c(:,:,cnt) = iffts(fft(iffts(A1b(:,:,cnt),[ 2]),[],2),[ 2]);
end;

%% keyboard;

for cnt=1:c,
  
  if pm == 1,
    Z1 = zeros(m,n);  Z1(1:2:end) = A1c(:,:,cnt);  
    Z2 = zeros(m,n);  Z2(2:2:end) = A2c(:,:,cnt); 
    Z1 = ffts(ifft2(ffts(Z1,[1 2])),[1 2]);
    Z2 = ffts(ifft2(ffts(Z2,[1 2])),[1 2]);  
    mask = abs(Z1)>0.1*max(abs(Z1(:)));
    indx=(mask~=0);

    ZZ = zeros(size(Z1));   
    ZZ( indx ) = angle(Z1(indx)./Z2(indx));
    phz = mean([ ZZ( ZZ>0); -ZZ(ZZ<0) ]);
    if (y(1,cnt)<0); phz=-phz; end;
  elseif ~isreal( y(2,cnt) ),
    phz = angle(y(2,cnt));
  else,
    phz = y(2,cnt);
  end;
  A1c(:,:,cnt) = A1c(:,:,cnt) .* exp(-j*phz/2); % angle(y(2,cnt)));
  A2c(:,:,cnt) = A2c(:,:,cnt) .* exp(j*phz/2); % angle(y(2,cnt)));
end;

a_out = A1c; % zeros(size(a));
% a_out(1:2:end,:,:) = A1c;
% a_out(2:2:end,:,:) = A2c;

b_out = A2c; % zeros(size(a));

