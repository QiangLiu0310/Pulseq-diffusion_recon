function lhdl = pltcmplx( in1, in2, in3, lbl )
%
% PLTCMPLX   sets up a GUI to compare the columns of complex matrices
%
% pltcmplx(a,b) or pltcmplx(a,b,c)
%
%

if ~ischar(in1),

  in1 = squeeze(in1); 
  fig = figure;
  if ~isreal(fig)                     % check if 'fig' is a matlab.ui.Figure class object
    fig = fig.Number;
  end;

  set(fig,'units','normalized');
  ans = get(fig,'position');
  set(fig,'position',[ 0.01 ans(2) 3*ans(3) ans(4) ]);
  figstr = num2str(fig);
  
  if nargin == 1, in2 = in1; end;
  in2 = squeeze(in2);
  Apltx = in1; Bpltx = in2;
  if nargin == 3, 
    if iscell(in3),
      lbl=in3; Cpltx = [];
    else,
      in3 = squeeze(in3); 
      Cpltx = in3; 
    end;
  else
    Cpltx = [];
  end;
  if nargin == 4, Cpltx = in3; end;
  
  [n,m] = size(Apltx);
  if ( (size(Bpltx,1) ~= n) | (size(Bpltx,2) ~= m) )
    close(fig);
    error('Input matrices must be the same size');
  end;
  
  [tmp,cnt] = max( max( abs(Apltx) ) );

  if m > 1,
    shdl = uicontrol('style','slider','position', [ 30 10 180 20 ], ...
                     'max', m, 'min', 1, 'value', cnt, ...
                     'tag', 'hdl1', ...
                     'tooltip','click here to advance through lines', ...
                     'callback', [ 'pltcmplx(''redraw'',' figstr, ',1)' ] );
    if m > 3,
      set(shdl, 'sliderstep', [ 1/m min(3/m,1) ] );
    end;
  else
    shdl = uicontrol('style','text','position', [ 30 10 180 20 ], ...
                     'value', 1, ...
                     'tag', 'hdl1' );
  end;
  
  vhdl = uicontrol( 'style', 'edit', 'string', num2str(cnt), ...
                    'position', [ 210 10 60 20 ], 'fontsize', 8, ...
		    'callback', [ 'pltcmplx(''text'',' figstr, ',2)' ], ...
		    'tooltip','Current line number', ...
                    'tag', 'hdl2' );
  dhdl = uicontrol( 'style','checkbox', 'string', 'fix axis', ...
                    'value',0, 'tag', 'hdl3', ...
                    'position', [ 280 10 80 20 ], ...
                    'tooltip', ...
                    'Set to plot all lines with same axis scale', ...
		    'callback', [ 'pltcmplx(''redraw'',' figstr, ',1)' ] );
  mhdl = uicontrol( 'style','checkbox', 'string', 'Mag/Phase', ...
                    'value',0, 'tag', 'hdl4', ...
                    'position', [ 370 10 80 20 ], ...
                    'tooltip', 'Plot Magnitude/Phase plots rather than Real/Imaginary data', ...
		    'callback', [ 'pltcmplx(''redraw'',' figstr, ',1)' ] );
  uhdl = uicontrol( 'style','checkbox', 'string', 'unwrap phase', ...
                    'value',0, 'tag', 'hdl5', ...
                    'position', [ 370 30 80 20 ], ...
                    'visible', 'off', ...
                    'tooltip', 'control unwrap/wrap phase plot', ...
		    'callback', [ 'pltcmplx(''redraw'',' figstr, ',1)' ] );

  tmp = '';

  for cnt = 1:(nargin-exist('lbl')),
    if exist('lbl'), 
      n_in = lbl{cnt};
    else,
      n_in = inputname(cnt);
    end;
    if strcmp( n_in, ''), n_in = sprintf('in%d',cnt); end;
    inlabels{cnt} = n_in;
    tmp = [ tmp n_in ];
    if cnt < nargin, tmp = [ tmp ' and ']; end;
  end;
  set(fig, 'Name', [ mfilename '(' tmp ')' ])
  
  set(fig, 'userdata', ...
	   struct( 'Apltx', Apltx, ...
		   'Bpltx', Bpltx, ...
		   'Cpltx', Cpltx, ...
		   'n', n, ...
		   'inlabels', { inlabels }, ...
		   'rmax', max([ real(Apltx(:)); real(Bpltx(:)) ]), ...
		   'rmin', min([ real(Apltx(:)); real(Bpltx(:)) ]), ...
		   'imax', max([ imag(Apltx(:)); imag(Bpltx(:)) ]), ...
		   'imin', min([ imag(Apltx(:)); imag(Bpltx(:)) ]), ...
		   'bwpr', 0, ... 
		   'shdl', shdl, ... 
		   'vhdl', vhdl, ... 
		   'dhdl', dhdl, ... 
		   'mhdl', mhdl, ... 
		   'uhdl', uhdl ) );
  
  % matlab's implementation of structured areas is soooo broken.  using
  % inlabels as a cell of strings causes the entire structure to be
  % copied to the largest dimension of inlabels.  What I want is Perl's
  % hash of hashes, where the data can be a pointer to more data !!!!!
  
else
  fig = in2;

  shdl = getfield( get(fig,'userdata'), 'shdl' );
  vhdl = getfield( get(fig,'userdata'), 'vhdl' );
  dhdl = getfield( get(fig,'userdata'), 'dhdl' );
  mhdl = getfield( get(fig,'userdata'), 'mhdl' );
  uhdl = getfield( get(fig,'userdata'), 'uhdl' );

  inlabels = getfield( get(fig,'userdata'), 'inlabels' ); 
  if ~iscell(inlabels), 
    tmp = {getfield( get(fig,'userdata'), 'inlabels' )}; 
    inlabels = tmp; 
  end;

end;

if strcmp( in1, 'text'),
  cnt = round( str2num( get(vhdl,'string') ) );
  if ~isempty( shdl )
    cnt( cnt < get(shdl,'min') ) = get(shdl,'min');
    cnt( cnt > get(shdl,'max') ) = get(shdl,'max');
    set( shdl, 'value', round(cnt) );
  else;
    cnt = 1;
  end;
end;
cnt = round( get( shdl, 'value') );
set( vhdl, 'string', num2str(round(cnt)) );

linenum = cnt;

rmax = getfield( get(fig,'userdata'), 'rmax' ); 
rmin = getfield( get(fig,'userdata'), 'rmin' ); 
imax = getfield( get(fig,'userdata'), 'imax' ); 
imin = getfield( get(fig,'userdata'), 'imin' ); 

n = getfield( get(fig,'userdata'), 'n' ); 
a = getfield( get(fig,'userdata'), 'Apltx' ); a = a(1:n,cnt); 
b = getfield( get(fig,'userdata'), 'Bpltx' ); b = b(1:n,cnt); 
if prod( size( getfield(get(fig,'userdata'), 'Cpltx') ) ),
  c = getfield( get(fig,'userdata'), 'Cpltx' ); c = c( 1:n,cnt );
else
  c = [];
end;

tmp = ''; lgnd_txt = '';
if ~isempty(inlabels)
 tmp = ' of ';
  for cnt = 1:length(inlabels),
    tmp = [ tmp inlabels{cnt} ]; 
    lgnd_txt = [ lgnd_txt '''' inlabels{cnt} ' ''' ];
    if cnt < length(inlabels), 
      tmp = [ tmp ' and ']; 
      lgnd_txt = [ lgnd_txt ', ' ]; 
    end;
  end;
end;

xpnts = (1:length(a))-length(a)/2;

if ( get( mhdl, 'value' ) )
  
  set( uhdl, 'visible', 'on' );

  dat1 = [ abs(a) abs(b) ];
  if ~isempty(c)
    dat1 = [ dat1 abs(c) ];
  end;
  ttl1 = [ 'Magnitude values' tmp ];
  
  if get(uhdl,'value'),
    dat2 = [ unwrap(angle(a)) unwrap(angle(b)) ]/pi;
    if ~isempty(c)
      dat2 = [ dat2 unwrap(angle(c))/pi ];
    end;
  else,
    dat2 = [ angle(a) angle(b) ]/pi;
    if ~isempty(c)
      dat2 = [ dat2  angle(c)/pi ];
    end;
  end;
  ttl2 = [ 'Phase values' tmp ]; 

  dat3 = [ angle(a./b) ]/pi;

  rmax = sqrt( max(rmax,abs(rmin))^2 + max(imax,abs(imin))^2 );
  rmin = 0;
  
  imax = max([ dat2(:);  1 ]);
  imin = min([ dat2(:); -1 ]);
  
else
  set( uhdl, 'visible', 'off' );
  dat1 = [ real(a) real(b) ];
  ttl1 = [ 'Real components' tmp ] ;
  dat2 = [ imag(a) imag(b) ];
  ttl2 = [ 'Imaginary components' tmp ];
  if ~isempty(c)
    dat1 = [ dat1 real(c) ];
    dat2 = [ dat2 imag(c) ];
  end;

  dat3 = [ angle(a./b) ]/pi;

end;  

if ~isreal( getfield( get(fig,'userdata'), 'Apltx') ),
  subplot(311);
end

if ( getfield( get(fig,'userdata'), {1}, 'bwpr')  > 0), 
  set(gca,'nextplot','replacechildren');
  set(gca,'colororder',[0 0 0]);
  set(gca,'linestyleorder', { 'x-'; 'o-.'; 'd-'});
end;
axs1 = axis;
plot( xpnts, dat1 );
if (sum( axs1 == [ 0 1 0 1]) == 4),
  axis('auto');
else;
  axs2 = axis;
  axis([ axs1(1:2) axs2(3:4) ]);
end;

if ( get(dhdl,'value') ),
  axis tight;
  %   axis([ 1 n rmin rmax ]);
end;
title(ttl1);
zoom off; 

% legend text set earlier.
% eval([ 'legend(' lgnd_txt ',0);' ]);
grid on;
for cnt=1:length(inlabels),
%  set( findobj('string', inlabels{cnt}), 'interp','none');
  % set( lhdl(cnt), 'interp', 'none');
end;

if ~isreal( getfield( get(fig,'userdata'), 'Apltx') ),
  subplot(312);

  if ( getfield( get(fig,'userdata'), {1}, 'bwpr')  > 0), 
    set(gca,'nextplot','replacechildren');
    set(gca,'colororder',[0 0 0]);
    set(gca,'linestyleorder', { 'x-'; 'o-.'; 'd:'});
  end;
  axs1 = axis;
  plot( xpnts, dat2 );
  if (sum( axs1 == [ 0 1 0 1]) == 4),
    axis('auto');
  else;
    axs2 = axis;
    axis([ axs1(1:2) axs2(3:4) ]);
  end;

  if ( get(dhdl,'value') ),
    axis tight; 
    % axis([ 1 n imin imax ]);
  end;
  title(ttl2);
  grid on;
  
  if ( get( mhdl, 'value' ) )
    ylabel( '\pi radians' );
  end;

  subplot(313);
  plot( xpnts, dat3 );

  if (sum( axs1 == [ 0 1 0 1]) == 4),
    axis('auto');
  else;
    axs3 = axis;
    axis([ axs1(1:2) axs3(3:4) ]);
  end;
  if ( get(dhdl,'value') ),
    axis tight;
    % axis([ 1 n imin imax ]);
  end;
  title(' Phase difference ');
  grid on;
  xlabel(['line ' num2str(linenum)]);
  
  if ( get( mhdl, 'value' ) )
    ylabel( '\pi radians' );
  end;

end;
