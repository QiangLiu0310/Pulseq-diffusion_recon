function a = fif(a)

  a = ffts(ifft2(ffts( a, 1:2)),1:2);
