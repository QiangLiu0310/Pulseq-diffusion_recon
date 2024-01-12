clc;close all;clear all;
A = zeros([96 96 4 4 ]);
A(1:4:end,:,:,:) = pro1_2(1:2:end,:,:,:,1);
A(2:4:end,:,:,:) = pro1_2(1:2:end,:,:,:,2);
A(3:4:end,:,:,:) = ifi(tmult(fif(pro3_4(2:2:end,:,:,:,1)),diag(exp(-j*4*pi*(1:96)/96)),2));
A(4:4:end,:,:,:) = ifi(tmult(fif(pro3_4(2:2:end,:,:,:,2)),diag(exp(-j*4*pi*(1:96)/96)),2));
imshow3(sos(fft2c(A),3))

B = zeros([96 96 4 4 ]);
B(1:4:end,:,:,:) = ifi(tmult(fif(flipdim(pro3_4(1:2:end,:,:,:,1),2)),diag(exp(-j*4*pi*(1:96)/96)),2));
B(2:4:end,:,:,:) = ifi(tmult(fif(flipdim(pro3_4(1:2:end,:,:,:,2),2)),diag(exp(-j*4*pi*(1:96)/96)),2));
B(3:4:end,:,:,:) = flipdim(pro1_2(2:2:end,:,:,:,1),2);
B(4:4:end,:,:,:) = flipdim(pro1_2(2:2:end,:,:,:,2),2);

