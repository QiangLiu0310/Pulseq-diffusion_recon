function a = fif3d(a)

  a = ffts(ifft(ifft2(ffts( a, 1:3)),[],3),1:3);
