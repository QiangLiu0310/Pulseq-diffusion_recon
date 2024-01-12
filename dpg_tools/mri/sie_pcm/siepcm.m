function [x,y,verbose] = siepcm(a,b,verbose,threshold,radius)
% SIEPCM - estimation of sub-pixel shifts between images using 
%          SVD of phase correlation 
%
% Inputs: 
%        a, b - input image data, in Fourier domain with DC centered
%     verbose - if non-zero, display debug plots and images
%    thrshold - a percentage of maximum value, used to set initial mask
%      radius - the width of data to use for LMS fit of phase slopes (<=0.25)

if nargin < 3, verbose = 0;     end;    % plot intermediate images
if nargin < 4, threshold = 0.0015; end;    % a percentage of max to set mask
if nargin < 5, radius = 0.25;   end;    % size of post-mask radius (<= 0.25)


[m,n] = size(a);

s = abs( a .* conj(b) ); s( s == 0 ) = 1;
c = (a .* conj(b)) ./ s; % abs( a .* conj(a) + 1e-10 );

s = ones(size(c));
s( abs(c) ~= 0 ) = abs(c( abs(c) ~= 0 ));


%%%% apply pre-masking here.

if exist('radius'), r=radius(1); else, r = 0.3;  end;

mask = zeros(size(c));
t = max(abs(a(:))).*threshold;

mask( abs(a) > t ) = 1;
mask = medfilt2( mask, [10 10] );       % 2D median filter, from image
                                        % processing toolbox
if (sum(mask(:)) == 0);
  error(sprintf([' mask in empty.  Need to lower the threshold (currently=' ...
  ' %f)'], threshold)); 
end;
                                        
[u d v]=svds(mask.*c,1);
% [u d v]=svds( c,1);

r2 = min(r,0.25);
% solve an LMS fit to the angles of the rank one matrix components...
if n > 3,
  n2 = length(v);
  %% apply post-masking:
  t_n = ceil( (0.5-r2)*length(v) ):floor( (0.5+r2)*length(v) );

  tmp = unwrap( angle(v(t_n)) ); A = [ (t_n-1)' ones(length(t_n),1) ];
  sys_y = A \ tmp;
  y = sys_y(1) * n / (2*pi);
  fullphz1 = unwrap(angle(v));
  fullphz1 = fullphz1 - fullphz1(floor(n2/2)) + ([floor(n2/2)-1 1]*sys_y);
  pltcmd1 = '1:length(v), fullphz1, t_n, [ unwrap( angle(v(t_n)) ) ( sys_y(1) * (t_n-1)'' + sys_y(2) ) ]';
else,
  pltcmd1 = ' 0,0 '; sys_y = [ 0 0]; y = 0;
end;

if m > 3,
  m2 = length(u);

  %% apply post-masking:
  t_m = ceil( (0.5-r2)*length(u) ):floor( (0.5+r2)*length(u) );
  
  tmp = unwrap( angle(u(t_m)) ); A = [ (t_m-1)' ones(length(t_m),1) ];
  sys_x = A \ tmp;
  x = sys_x(1) * m / (2*pi);
  fullphz2 = unwrap(angle(u));
  fullphz2 = fullphz2 - fullphz2(floor(m2/2)) + ([floor(m2/2)-1 1]*sys_x);
  pltcmd2 = '1:length(u), fullphz2, t_m, [ unwrap( angle(u(t_m)) ) ( sys_x(1) * (t_m-1)'' + sys_x(2) ) ]';
else,
  pltcmd2 = ' 0, 0 '; sys_x = [ 0 0]; x = 0;
end;

% if verbose mode is on, show the estimation quality images
if verbose,  
  figure(1);
  disp([ 'plot( ' pltcmd1 ', ''o-'', ' pltcmd2 ', ''*-'');' ]);
  eval([ 'plot( ' pltcmd1 ', ''o-'', ' pltcmd2 ', ''*-'');' ]);
  title(' plot of least squares linear fit ');
  
  
  Q = exp( j*2*pi* x * [0:m-1]'/m ) * exp( -j*2*pi* y * [0:n-1]/n );
  q = u*v'; 
  
  % phase align each of these beasts so that the abs( c - Q ) looks good.
  mc = floor(m/2); nc = floor(n/2);
  Q = Q * exp(angle(Q(mc,nc)'*c(mc,nc))*j);
  q = q * exp(angle(q(mc,nc)'*c(mc,nc))*j);
  figure(2);
  imagesc( angle([ c mask.*c; u*v' Q ]) ); colormap('hsv');
  title('[ Q | mask \circ Q ; u v^H | est of (u v^H) ]');
  
  disp([ '    x_0: ' num2str(x) '       y_0: ' num2str(y) ]);
  keyboard;
end;

% Copyright 2002-2005 William Scott Hoge (shoge at ieee dot org) 
% All rights reserved. Licensed according to the GNU General Public Licence
% (see http://www.gnu.org/copyleft/gpl.html)