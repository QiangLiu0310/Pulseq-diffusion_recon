function res = fft3call(x)
fctr = size(x,1)*size(x,2)*size(x,3);

X = fftshift(fft(ifftshift(x,1),[],1),1);

X = fftshift(fft(ifftshift(X,2),[],2),2);

res = fftshift(fft(ifftshift(X,3),[],3),3) / sqrt(fctr);



