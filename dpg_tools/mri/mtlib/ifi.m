function a = ifi(a)
  a = iffts(fft2(iffts( a, 1:2)),1:2);
