function [k_new,Is] = reconhd(pksp,Np,Nm)
%
% port of pifft.m to non-power-of-2 zero-padded data
%

% Np  - number of phase encodes
% Nm  - full span of image to reconstruct
% Npp - number of points in conjugate-symmetric region 

sz = size(pksp);

if nargin==1
  Nm = size(pksp,1); % /2;
  Np = (max(find(pksp(:,1))));
end

Npp = 2*( (Np)-Nm/2 );



% Convert data into hybrid space by inverse FFT along FE dir.
% for better numerical accuracy.
pksp = ifft(pksp,[],2);

% Frequency Encode direction is assumed to be the longer direction
[tmp Nf] = size(pksp);

% Number of phase encodes in the full k-space
% Np = 2^nextpow2(Npp+1);

% Number of high frequency phase encodes
NH = Np - Npp;

% NUmber of low frequency phase encodes
NL = Np - NH;

% Row indices of high and low frequency lines
HFind = 1:NH;
LFind = (NH+1):(NH+NL);

% Create weighting window
%
% It multiplies all un-partnered high frequency lines by 2.0
% Low frequency lines are multiplied by a ramp that approaches 
% zero at the zero-padded edge of k-space.  The weighting
% factors of two partnered low frequency lines should have a 
% sum of 2.0
w = zeros(Nm,1);
w(HFind) = 2;
rstep = 2/(NL+1);
w(LFind) = [(2-rstep):-rstep:rstep];

% Create weighted partial k-space
HFksp = zeros(Nm,Nf);
HFksp(1:Np,:)=pksp(1:Np,:);
HFksp = HFksp.*repmat(w,1,Nf);

% Create Low Frequency k-space 
LFksp = zeros(Nm,Nf);
LFksp(LFind,:) = pksp(LFind,:);

% Low frequency image
Rc = ifft(LFksp,[],1);

% Unsynchronous image
Ic = ifft(HFksp,[],1);


% Synchronous image
Is = Ic.*exp(-i*angle(Rc));


k_new = fftshift(fft2(Is));

% Demodulated image
% m = real(Is);

