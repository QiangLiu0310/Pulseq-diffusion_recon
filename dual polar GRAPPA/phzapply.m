function [a_out] = phzapply(a,y);
% [a_out] = phzapply(a,y);
%
% apply the EPI phase shift parameters calcuated by phzshift or
% comp_local_pc to the original k-space data
%


[m n c]=size(a);

[pm pc ] = size(y);

A1 = a(1:2:end,:,1:c);
A1a = ffts(ifft(ffts(A1,[ 2]),[],2),[ 2]);

A2 = a(2:2:end,:,1:c);
A2a = ffts(ifft(ffts(A2,[ 2]),[],2),[ 2]);

A1b = zeros(size(A1a));
A2b = zeros(size(A2a));
for cnt=1:c,
  A1b(:,:,cnt) = A1a(:,:,cnt) * diag( exp( +j*y(1,cnt)/2*((1:size(A2a,2))-size(A2a,2)/2)));
  A2b(:,:,cnt) = A2a(:,:,cnt) * diag( exp( -j*y(1,cnt)/2*((1:size(A2a,2))-size(A2a,2)/2)));
end;

A1c = iffts(fft(iffts(A1b,[ 2]),[],2),[ 2]);
A2c = iffts(fft(iffts(A2b,[ 2]),[],2),[ 2]);

%% keyboard;

for cnt=1:c,
  if ~isreal( y(2,cnt) ),
    phz = angle(y(2,cnt));
  else,
    phz = y(2,cnt);
  end;
  A1c(:,:,cnt) = A1c(:,:,cnt) .* exp(-j*phz/2); % angle(y(2,cnt)));
  A2c(:,:,cnt) = A2c(:,:,cnt) .* exp(j*phz/2); % angle(y(2,cnt)));
end;

a_out = zeros(size(a));
a_out(1:2:end,:,:) = A1c;
a_out(2:2:end,:,:) = A2c;

