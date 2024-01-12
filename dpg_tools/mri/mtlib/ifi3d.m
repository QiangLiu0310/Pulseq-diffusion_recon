function a = ifi3d(a)

  a = iffts(fft(fft2(iffts( a, 1:3)),[],3),1:3);
