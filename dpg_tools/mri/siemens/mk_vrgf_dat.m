function vrgf = mk_vrgf_dat( inres, outres, t_ru, t_sd, t_ft, t_sw )
%
% vrgf = mk_vrgf_dat( inres, outres, t_ru, t_sd, t_ft, t_sw )
%
% k : input k-space
%
% freqres : output resolution in the frequency dimension
% t_ru : time, ramp up (attack)  [ms]
% t_sd : time, sampling delay    [ms]
% t_ft : time, flat-top          [ms]
% t_sw : sampling window         [ms]
%
% we assume that the pulse is symmeteric, e.g. ramp down time, t_rd = t_ru
%

t_rd = t_ru;

%
% total area of ramp: 1/2 * Amp * t_ru
%

% plot out the trapezoidal gradient pulse, in msec
grad = [ ( 0 : (t_ru) ) ...
         t_ru*ones(1,t_ft-1) ...
         (t_ru-[0:t_rd]) ];
gradpts = [ 0 : length(grad)-1 ] - (length(grad)-1)/2;


% use the 'sampling delay' parameter to align the sampling window to the gradient pulse, 
sampwin = zeros(1,length(grad));
sampwin( t_sd+(0:t_sw) ) = 1;



% calculate the sampling time point values,  sampling the acquisition uniformly.
% tpts    = t_sd + ( (t_sw/(inres-1))/2 : t_sw/(inres+1) : t_sw-(t_sw/(inres-1))/2 );
samplewidth = t_sw/(inres);
tpts  = -t_sw/2+samplewidth/2 : samplewidth : t_sw/2-samplewidth/2;

% calculate the regions for each part of the gradient
indx1 = find( tpts < -t_ft/2 );                  % 1) ramp up
indx2 = find( (tpts>=-t_ft/2)&(tpts<=t_ft/2) );  % 2) flat-top
indx3 = find( tpts > t_ft/2 );                   % 3) ramp down

% for each region, determine the gradient value at each sampling time point
gdpts = zeros(1,inres); 
gdpts(1:length(indx1)) = t_ru - ( -t_ft/2 - tpts(indx1) );  
gdpts(length(indx1)+[1:length(indx2)]) = t_ru; 
gdpts(length(indx1)+length(indx2)+[1:length(indx3)]) = t_ru - [tpts(indx3) - t_ft/2];

%%% plot the result.
% plot( gradpts, grad, tpts, gdpts, 'kx' );
% newpts  = 0 : t_sw/(outres-1) : t_sw;

%%% calculate the relative gradient applitude change at each sampled time point
z  = cumsum( gdpts ) - gdpts(1);
z2 = [ z*(inres-1)/max(z(:)) - inres/2 + 0.5 ];   % ... this gives the k-space location at each time point

Z = zeros([inres outres]);
W = zeros([inres outres]);

stp = (inres-1)/(outres-1);
for cnt=1:outres,                       % for each output point ...
  
  % form a sinc interpolator:
  Z(:,cnt) = sinc( (z2 - ((cnt-1)*stp - (inres-1)/2) ) * 0.6 )'; % 76/inres )' ; 
  % and a blackman window, centered at the sampling point
  % tmp = ( 0.42 - 0.5*cos(2*pi*(z2 - ((cnt-1)*stp - (inres-1)/2) )) + 0.08*cos(4*pi*(z2- ((cnt-1)*stp - (inres-1)/2) )) );
  tmp = ( 0.42 - 0.5*cos(2*pi*(z2-(cnt-1)*inres/outres)/(inres-1)) + 0.08*cos(4*pi*(z2-(cnt-1)*inres/outres)/(inres-1)) );

  % the blackman window is cyclic, so clear out the values that appear beyond
  % the primary bump
  indx = floor(inres/2+(cnt*inres/outres));
  tmp( indx:end ) = 0;
  indx2 = 1:(indx-inres);
  tmp( indx2(find(indx2)>0) ) = 0;

  W(:,cnt) = tmp(:);
end;
  
% form the resampling matrix, by applying the windowing function to the sinc interpolator
vrgf = W.*Z;

% scale the coefficients, so that the max DC component is == 1
vrgf = vrgf * diag( 1./sum(vrgf,1) );


%%%
% Note:  The filter that is generated here is not quite correct.  It should be symmetric, s.t.
% vrgf == rot90(rot90(vrgf));
%
% to correct, compute the phase shift between the middle two columns:

p = comp_local_pc( fif(vrgf(:,outres/2)), fif(flipdim(vrgf(:,1+outres/2),1)) );

% and then apply 1/2 of the shift to the resampling operator
vrgf = exp(-j*p(2)/2) * ...
      iffts(fft(iffts( ...
          tmult( ffts(ifft(ffts(vrgf,1),[],1),1), diag(exp( +j*p(1)/2*( (1:size(vrgf,1))-size(vrgf,1)/2) )), 1 ), ...
	    1),[],1),1);

%% p.s. there is probably a way to do this above... e.g. by changing
%% one of the step parameters into a partial step.
