function result = im(in1,fig,parm1,parm2)
% IM  show image(s) of complex valued MRI data
%
% usage: im( A, [fig, settings ])
%
%     A         -  image matrix to show (supports upto 4D matrices)
% 
%     optional parameters:
%     fig       -  either the figure number of a figure label.
%                  if a label is given, the image will appear in a
%                  currently open figure with that label.  If no such
%                  figure exists, then a new figure is created and given
%                  that label.
%
%     settings  -  initially show the image by one of the following
%                  types: {'magnitude','phase','real','imag','gop'}
%

if nargin < 3, parm1 = []; end;
if nargin < 4, parm2 = []; end;

if ~ischar(in1),
  
  if ( nargin < 2 ); fig = []; end;

  tagstr = [];
  if isempty(fig),
    fig = figure;
    if ~isreal(fig) % check if 'fig' is a matlab.ui.Figure class object
      fig = fig.Number;
    end;
  elseif ischar(fig),
    tagstr = fig;
    fig = findobj('tag',tagstr);
    if isempty(fig);
      fig = figure;
      set(fig,'tag',tagstr);
    end;
  else
    if ~isreal(fig)                     % check if 'fig' is a matlab.ui.Figure class object
      fig = fig.Number;
    end;
    figure(fig)
  end;
  
  zoom off;
  c = [];
  usrdat = get(fig,'userdata');
  if ( ~isempty( usrdat ) ),
    if ( isfield( usrdat, 'cmap' ) )
      c = getfield( get(fig,'userdata'), 'cmap' );
    else
      c = [];
    end;
    clf;
%     hdls = { 'mhdl', 'thdl', 'ehdl', 'chdl', 'shdl', 'vhdl', 'ghdl' };
%     for length(
%     if ~isempty( findobj(get(fig,'children'), 'tag','hdl4') )
%       delete( getfield( get(fig,'userdata'), 'mhdl' ) );
%     end;

%     delete( getfield( get(fig,'userdata'), 'thdl' ) );
%     delete( getfield( get(fig,'userdata'), 'ehdl' ) );
%     delete( getfield( get(fig,'userdata'), 'shdl' ) );
%     delete( getfield( get(fig,'userdata'), 'vhdl' ) );
%     delete( getfield( get(fig,'userdata'), 'chdl' ) );
%     delete( getfield( get(fig,'userdata'), 'ghdl' ) );
  end;

  msk = [];
  if ~isempty(parm2)
    if ( sum(size(in1) == size(parm2)) == length(size(in1)) )
      msk = parm2;
    end;
  end;
    
  cb = [ mfilename '(''redraw'', ', num2str(fig), ')' ];

  %%% Wierd flickering funniness with the 'popup' style icons.  
  %%% so, set the 'invisible' ...
    mhdl = uicontrol( 'style','popupmenu', ...
                      'visible', 'on', ...
                      'string', {'magnitude','phase','real', 'imag', 'gop' }, ...
		      'value',1, 'tag', 'mhdl', ...
		      'position', [ 10 10 80 20 ], ...
		      'tooltip', ...
		      'Switch display between Magnitude, Phase, Real, Imaginary, and GOP images', ...
		      'callback', cb );
    if ischar( parm1 ), 
      tmp = find( strcmp( get( mhdl, 'string' ), parm1 ) );
      if tmp ~= 0,
        set(mhdl,'value',tmp);
      end;
    end;

    %%% and use Menu settings instead to update the value field(s)
    mhdl2 = [];
%    mhdl2 = addmenu(fig,7, ...
%                    { '&Display mode', '&magnitude', '&phase', '&real', '&imag', '&gop'}, ...
%                    { '', ...
%                     [ 'im(''setdisp'',' num2str(fig) ',1)'] ...
%                     [ 'im(''setdisp'',' num2str(fig) ',2)'] ...
%                     [ 'im(''setdisp'',' num2str(fig) ',3)'] ...
%                     [ 'im(''setdisp'',' num2str(fig) ',4)'] ...
%                     [ 'im(''setdisp'',' num2str(fig) ',5)'] ...
%                    } );
%    set(mhdl2(2),'checked','on');

  
  
  thdl = uicontrol( 'style','checkbox', 'string', 'Threshold', ...
                    'value',0, 'tag', 'thdl', ...
                    'position', [ 10 30 80 20 ], ...
                    'tooltip', 'Turn thresholding on or off', ...
                    'callback', cb );

  ehdl = uicontrol( 'style','checkbox', 'string', 'Difference', ...
                    'value',0, 'tag', 'ehdl', ...
                    'position', [ 10 50 80 20 ], ...
                    'tooltip', 'compute difference between left half and right half', ...
                    'callback', cb );
% 
%   ehdl = uicontrol( 'style','checkbox', 'string', 'Overlay', ...
%                     'value',0, 'tag', 'ohdl', ...
%                     'position', [ 10 70 80 20 ], ...
%                     'tooltip', 'overlay the left half of the image (in color) ontop of the right half (in gray)', ...
%                     'callback', cb );

  chdl = uicontrol( 'style','checkbox', 'string', 'common caxis', ...
                    'value',0, 'tag', 'chdl', ...
                    'position', [ 10 90 80 20 ], ...
                    'tooltip', 'toggle caxis from ''auto'' to max(in1)', ...
                    'callback', cb );

  in1 = squeeze(in1);
  if (length(size(in1)) > 4),
    error([ '''im'' can only handle upto 4D data. '...
              'use ''squeeze'' to reduce the number of dimensions.' ]);
  end;
  
  shdl = []; vhdl = []; mvhdl = [];
  if length( size( in1 ) ) >= 3,
    shdl = uicontrol('style','slider','position', [ 90 20 180 20 ], ...
                     'max', size(in1,3), 'min', 1, 'value', 1, ...
                     'tag', 'shdl', ...
                     'tooltip','click here to advance through index 3', ...
                     'sliderstep', [ 1 10 ]/(size(in1,3)-1), ...
		     'callback', cb );
    vhdl = uicontrol( 'style', 'text', 'string', num2str(1), ...
                      'position', [ 270 20 60 20 ], 'fontsize', 8, ...
                      'tooltip','Current line number', ...
                      'tag', 'vhdl' );
  
    mvhdl = uicontrol( 'style','checkbox', 'string', 'movie', ...
                       'tag','mvhdl', ...
                       'tooltip', 'play a movie of all frame', ...
                       'position', [ 10 170 80 20 ], ...
                       'callback', [ mfilename '(''playmovie'',' num2str(fig) ');' ] ...
                      );

  end;
  if length( size( in1 ) ) == 4,
    shdl(2) = uicontrol('style','slider','position', [ 90 0 180 20 ], ...
                        'max', size(in1,4), 'min', 1, 'value', 1, ...
                        'tag', 'shdl', ...
                        'tooltip','click here to advance through index 4', ...
                        'sliderstep', [ 1 10 ]/(size(in1,4)-1), ...
                        'callback', cb );
    vhdl(2) = uicontrol( 'style', 'text', 'string', num2str(1), ...
                         'position', [ 270 0 60 20 ], 'fontsize', 8, ...
                         'tooltip','Current line number', ...
                         'tag', 'vhdl' );
  end;

  rsoshdl = uicontrol( 'style', 'checkbox', 'string', 'RSoS', ...
                       'position', [340 0 100 20], ...
                       'tag','rsos','visible','off','callback', cb );
  if length( size( in1 ) ) > 2,
    set(rsoshdl,'visible','on');
  end;

  uicontrol( 'style', 'text', 'string', 'Display Range:', ...
             'position', [ 340 20 100 20 ], ... % 'fontsize', 8, ...
             'tag', 'drhdl' );
  rlhdl = uicontrol( 'style', 'edit', 'string', 0, ...
                     'position', [ 450 20 30 20 ], ... % 'fontsize', 8, ...
                     'tag', 'rlhdl', 'callback', cb );
  rhhdl = uicontrol( 'style', 'edit', 'string', 1, ...
                     'position', [ 490 20 30 20 ], ... % 'fontsize', 8, ...
                     'tag', 'rhhdl', 'callback', cb );


  
  if (0), % ~isempty(usrdat)
    mxt = getfield( get(fig,'userdata'), 'mxt' );
    mnt = getfield( get(fig,'userdata'), 'mnt' );
    mxi = getfield( get(fig,'userdata'), 'mxi' );
    mni = getfield( get(fig,'userdata'), 'mni' );
    mxp = getfield( get(fig,'userdata'), 'mxp' );
    mnp = getfield( get(fig,'userdata'), 'mnp' );
    mxr = getfield( get(fig,'userdata'), 'mxr' );
    mnr = getfield( get(fig,'userdata'), 'mnr' );
  else
    mxt = max( abs(double(in1(:)))); mxp = max(angle(double(in1(:))));
    mnt = min( abs(double(in1(:)))); mnp = min(angle(double(in1(:))));
    mxi = max(imag(double(in1(:)))); mxr = max( real(double(in1(:))));  
    mni = min(imag(double(in1(:)))); mnr = min( real(double(in1(:))));  
  end;
  
  
  ahdl = uicontrol( 'style','popupmenu', ...
                    'string', {'image', 'equal', 'normal'}, ...
                    'visible','off', ...
                    'tag', 'ahdl', ...
		    'position', [ 10 110 80 20 ] );

  cb2 = [ '''''toggle'''', ', num2str(fig)  ];
  cb3 = [ 'eval([ ''' mfilename '( ' cb2 ','' num2str(t) '')'' ])' ];

  
  snrhdl = uicontrol( 'style','pushbutton', 'string', 'calc SNR', ...
                      'position', [ 10 150 80 20 ], ...
                      'callback', [ mfilename '(''calcsnr'',' num2str(fig) ');' ] ...
                      );

  pnthdl = uicontrol( 'style','pushbutton', 'string', 'print all', ...
                      'position', [ 10 170 80 20 ], ...
                      'callback', [ mfilename '(''printall'',' num2str(fig) ');' ] ...
                      );


  
  ghdl = uicontrol( 'style','pushbutton', 'string', 'get point', ...
		    'position', [ 10 110 80 20 ], ...
		    'callback', ...
		    [ 'tmp=ginput(1); t = ( round(tmp(1)-1)*' ...
		      'getfield( get(gcf,''userdata''), ''size_X'',{1}) '...
		      '+ round(tmp(2)) ), ' cb3 ]); 
		      
  % a(t)= not(a(t));b(t)=not(b(t));' cb ] ); 
  % disp( [ floor(tmp(1)) getfield( get(gcf,''userdata''), ''size_X'',{1})' ...
  %	    ' floor(tmp(2)) ]) ; 
  
  if isempty(c),
    c = colormap;                         % fetch the current colormap
  else
  end;
  set( fig, 'userdata', struct( 'im_X', double(in1), ...
				'size_X', size(in1), ...
                                'mask', msk, ...
                                'size_mask', size(msk), ...
                                'cmap', c, ...
				'mxt', mxt, ...
				'mnt', mnt, ...
				'mxp', mxp, ...
				'mnp', mnp, ...
				'mxr', mxr, ...
				'mnr', mnr, ...
				'mxi', mxi, ...
				'mni', mni, ...
				'ahdl', ahdl, ...
				'ghdl', ghdl, ...
				'mhdl', mhdl, ...
				'menu', mhdl2, ...
				'thdl', thdl, ...
				'ehdl', ehdl, ...
				'vhdl', vhdl, ...
				'chdl', chdl, ...
				'rsoshdl', rsoshdl, ...
				'rlhdl', rlhdl, ...
				'rhhdl', rhhdl, ...
				'mvhdl', mvhdl, ...
				'shdl', shdl ) );

  if isempty(tagstr),
    set(fig, 'Name', [ mfilename '( ' inputname(1) ' )' ]);
  else,
    set(fig, 'Name', [ tagstr ] );
  end;
  in1 = 'redraw';
else,
  if ( nargin < 2 ); fig = gcf; end;

end;

mhdl = getfield( get(fig,'userdata'), 'mhdl' );
thdl = getfield( get(fig,'userdata'), 'thdl' );
ehdl = getfield( get(fig,'userdata'), 'ehdl' );
shdl = getfield( get(fig,'userdata'), 'shdl' );
vhdl = getfield( get(fig,'userdata'), 'vhdl' );
chdl = getfield( get(fig,'userdata'), 'chdl' );
ahdl = getfield( get(fig,'userdata'), 'ahdl' );
rsoshdl = getfield( get(fig,'userdata'), 'rsoshdl' );
rlhdl = getfield( get(fig,'userdata'), 'rlhdl' );
rhhdl = getfield( get(fig,'userdata'), 'rhhdl' );

% tmp = getfield( get(fig,'userdata'), 'im_X' ); 
%, {:,:,round(get(shdl,'value')} );

if ~isempty(shdl),
  tmp2 = [];
  % specify which rows and colums if im_X to read:
  cmd = [ ' im_X = getfield( get(fig,''userdata''), ''im_X'', ' ...
        ' { 1:getfield( get(fig,''userdata''), ''size_X'',{1}), ' ...
        ' 1:getfield( get(fig,''userdata''), ''size_X'',{2}), ' ];
  % specify which elements from the 3rd (and 4th) dimension to read:
  nar = length(vhdl) - get(rsoshdl,'value');

  for tmp=1:nar, % length(vhdl)
    tmp2 = round(get(shdl(tmp),'value'));
    set(vhdl(tmp),'string',num2str( round(get(shdl(tmp),'value')) ));
    cmd = [ cmd num2str(tmp2) ];
    if tmp<nar, cmd = [ cmd ', ' ]; end;
  end;
  if ( get(rsoshdl,'value') )
    cmd = [ cmd ' '':'' ' ];
    set( vhdl(end), 'visible', 'off');
    set( shdl(end), 'visible', 'off');
  else
    set( vhdl(end), 'visible', 'on');
    set( shdl(end), 'visible', 'on');
  end;
  cmd = [cmd ' } );' ];
  eval(cmd);                            % fetch the image data to display

  if ( sum( getfield( get(fig,'userdata'), 'size_mask' ) ) ~= 0 )
    tmp2 = [];
    cmd = [ ' msk = getfield( get(fig,''userdata''), ''mask'', ' ...
            ' { 1:getfield( get(fig,''userdata''), ''size_X'',{1}), ' ...
            ' 1:getfield( get(fig,''userdata''), ''size_X'',{2}), ' ];
    for tmp=1:length(vhdl)
      tmp2 = round(get(shdl(tmp),'value'));
      set(vhdl(tmp),'string',num2str( round(get(shdl(tmp),'value')) ));
      cmd = [ cmd num2str(tmp2) ];
      if tmp<length(vhdl), cmd = [ cmd ', ' ]; end;
    end;
    cmd = [cmd ' } );' ];
    eval(cmd);
  else
    msk = [];
  end;
else
  im_X = getfield( get(fig,'userdata'), 'im_X' );
  if ( sum( getfield( get(fig,'userdata'), 'size_mask' ) ) ~= 0 )
    msk = getfield( get(fig,'userdata'), 'msk' );
  else
    msk = [];
  end;
end;  

if ( get(rsoshdl,'value') )
  if ( get( mhdl, 'value' ) == 2 ), % complex data shown
    im_X = vbc( im_X );
  else
    im_X = sqrt(sum(abs( im_X ).^2,length(size(im_X))));
  end;
end;

if strcmp(in1,'redraw'),
  if isempty(mhdl),
    rsp = 1; 
  else,
    rsp = get( mhdl, 'value' );
  end;

  tmp2 = getappdata(fig,'ZoomOnState');
  if ( ~isempty( tmp2 ) ),
    cura = axis;
  end;

  % from CF Westin's gopimage.m
  magf = abs(im_X); % or use: exp(gamma*log(abs(f)+eps));
  maxf = max(abs( getfield( get(fig,'userdata'), 'im_X' ) ));
  maxf = max(maxf(:));
  mbits = 4; abits = 8-mbits;
    
  switch rsp,
   case 1,
    tmp = abs( im_X );  ttl = 'Intensity values';
    maxca = 'mxt'; minca = 'mnt';
   case 2,
    tmp = angle( im_X )./pi; ttl = 'Phase values';
    maxca = 'mxp'; minca = 'mnp';
   case 3,
    tmp = real( im_X );  ttl = 'Real values';
    maxca = 'mxr'; minca = 'mnr';
   case 4,
    tmp = imag( im_X );  ttl = 'Imaginary values';
    maxca = 'mxi'; minca = 'mni';
   case 5,
    % from CF Westin's gopimage.m
    % if (maxf == 0),
    %   fprintf('Zero magnitude. Image is rendered with max magnitude\n');
    %   magf = (2^mbits - 1) * ones(magf);
    % else
      magf = floor(magf / maxf * (2^mbits - 0.001));
    % end
    % Scale the argument to the interval 0 ... (2^abits - 1)
    argf = floor((angle(im_X * exp(-i * pi)) + pi) / (2 * pi) * (2^abits - 0.001));
    tmp = (magf * (2^abits) + argf + 1);
    ttl = 'CF''s gop style';
   
   otherwise,
    tmp = abs( im_X );  ttl = 'Intensity values';
  end;

  if ( get( ehdl, 'value' ) ),          % show difference between L and R
    t = size(tmp,2)/2;
    if (rsp == 2),
      tmp = angle( exp(j*pi*tmp(:,1:t)) ./ exp(j*pi*tmp(:,t+(1:t))) )/pi;
    else
      tmp = abs(tmp(:,1:t) - tmp(:,t+(1:t)));
    end;
    if ~isempty(msk)
      sz = size(msk);
      msk = max( reshape([ msk(:,1:t) msk(:,t+(1:t)) ],[sz(1) t 2]), [], 3 );
    end;
  else,
    t = size(tmp,2);
  end;
  
  if ( get( thdl, 'value' ) ),          
    % clip low values if threshold button is set to 'on'
    tmp( abs(im_X(:,1:t)) < 0.05*max(abs(im_X(:))) ) = min(tmp(:)) ;
  end;

  d = (max(tmp(:)) - min(tmp(:)));
  lo = str2num( get( rlhdl, 'string' )) * d + min(tmp(:));
  hi = str2num( get( rhhdl, 'string' )) * d + min(tmp(:));

  tmp( tmp > hi ) = hi;
  tmp( tmp < lo ) = lo;
  
  pppp = findobj(fig,'tag','main image');
  qqqq = findobj(fig,'tag','image mask');
  if isempty(pppp),
    pppp = image(tmp); % imagesc( tmp );
    set(pppp,'tag','main image');
  else,
    set(pppp,'cdata',tmp);
  end;
  if ( isempty(qqqq) && ~isempty(msk) ),
    hold on;
    qqqq = image(zeros([ size(tmp) 3]));
    set(qqqq,'tag','image mask');
    hold off;
    set(qqqq,'alphadata',msk);
  else
    set(qqqq,'cdata',zeros([ size(tmp) 3]));
    set(qqqq,'alphadata',msk);
  end;
  
  str = get(ahdl,'string');
  axis( str{ get(ahdl,'value') } );

  if ( ~isempty( tmp2 ) ),
    zoom('on');
    axis(cura);
%    keyboard;
  end;

  % if ( sum( tmp2 == [ 0 1 0 1 ] ) ~= 4 ),
  %   axis(tmp2);
  % end;
  % set( get(pppp,'parent'), 'XLim', tmp2(1:2), 'YLim', tmp2(3:4) )

  if ( get( chdl, 'value' ) & ( rsp ~= 5 ) ),
    cmn = getfield( get(fig,'userdata'), minca );
    cmx = getfield( get(fig,'userdata'), maxca );
    % disp(['using common color axis:' num2str([ cmn cmx ])]);
    caxis([ cmn cmx ]);
  else
    caxis('auto');
  end;
  
  curclim = [];
  if strcmp( get(gca,'climmode'), 'manual');
    curclim = get(gca,'clim');
  end;
  

  argstep = 256 / (2^abits);       % Step length to be used in goptab
  gtab    = goptab; % hsv(256);                % Load goptab colour table
  gtab    = [gtab;gtab(1,:)]; % Make it cyclic to allow interpolation for all angles
  argtab  = gtab(1:argstep:256,:); % Compile an argument colour table
  
  magstep = 255 / (2^mbits - 1);          % Step length to be used in gray
  gtab    = gray(256);                    % Load gray colour table
  magtab  = gtab(round(1:magstep:256),:); % Compile a magnitude colour table
  
  % Make a composite magnitude/colour table. Argument are the LSB
  C = [kron(magtab(:, 1), argtab(:, 1)) kron(magtab(:, 2), argtab(:, 2)) ...
       kron(magtab(:, 3), argtab(:, 3))];
  
  curcmap = colormap;
  if ( rsp == 5 ),                 % if using GOP, do things a bit different
    set(pppp,'CDataMapping','direct');
    if ( prod( double(size(curcmap) == size(C) ) ) ),
      if (sum( curcmap(:) - C(:) ) ~= 0),
        % disp('saving current colormap');
        img_data = get(fig,'userdata');
        img_data.cmap = curcmap;
        set( fig, 'userdata', img_data ); 
        % only save the previous colormap if it is not the gopmap
      end;
    else
      % disp('saving current colormap');
      img_data = get(fig,'userdata');
      img_data.cmap = curcmap;
      set( fig, 'userdata', img_data ); 
    end;
    colormap(C);

    % % construct a relevant colorbar
    % x = linspace(-1,1,35);
    % [X,Y]=ndgrid(x,x);
    % f = Y-j*X;
    % 
    % a = colorbar;
    % ans = get(a,'Position')
    % set(a,'position', [ ans(1) ans(2)+0.4 2*ans([3 3]) ]);
  else
    set(pppp,'CDataMapping','scaled');
    if ( prod( double(size(curcmap) == size(C) ) ) ),
      if (sum( curcmap(:) - C(:) ) == 0),
        % disp('current cmap == gop :: retrieve old cmap');
        colormap( getfield( get(fig,'userdata'), 'cmap' ) );
      end;
    end;
    % colorbar; 
  end;
  axis( get( findobj(fig,'tag','main image'), 'parent') );
  title(ttl);

  if (~isempty(curclim)),
    set(gca,'clim',curclim);
  end;

elseif strcmp( in1, 'calcsnr'),

  % fig = get(gcbo,'parent');

  %%% clean up the overlay graphics if still around...
  l1 = findobj(fig,'tag','box1'); if ~isempty(l1), delete(l1); end;
  l2 = findobj(fig,'tag','box2'); if ~isempty(l2), delete(l2); end;
  
  ihdl = findobj(fig,'tag','overlay');
  if ~isempty(ihdl), set(ihdl,'visible','off'); end;
  
  ia = findobj(fig,'tag','main image');
  la = get( get( ia, 'parent' ), 'xlabel' );
  set(la,'string','select a region with reasonably constant signal');


  %%% select region parameters
  t = -3:3;
  t1 = ones(length(t),1)*ginput(1)+[ t' t'];

  A = get(ia,'cdata'); % field( get(fig,'userdata'), 'im_X' );
  ovrlay = zeros(size(A));
  ovrlay( round(t1(:,2)), round(t1(:,1)) ) = 1;
  
  ihdl = findobj(fig,'tag','overlay');
  if isempty(ihdl);
    tmp = axes;
    ihdl = imagesc(ovrlay); axis('image');
    ahdl = findobj(fig,'tag','ahdl');
    str = get(ahdl,'string');
    axis( str{ get(ahdl,'value') } );

    set(get(ihdl,'parent'),'visible','off');
  end;
  set(ihdl,'cdata',ovrlay,'alphadata',0.7*ovrlay,'tag','overlay','visible','on');
  
  cnr = vec( round([ t1(1,:) ; t1(size(t1,1),:) ]) + [ -0.5 -0.5; 0.5 0.5 ] );
  
  axis(get(ihdl,'parent'));
  l1 = line([ cnr(1) cnr(1) cnr(2) cnr(2) cnr(1) ], ... 
            [ cnr(3) cnr(4) cnr(4) cnr(3) cnr(3) ], ... 
            'color',[ 1 1 1 ], 'tag', 'box1');
  
  
  % A2 = ifft2(fftshift(fft2(A2))
  % [ sqrt( A2(:)'*A2(:) )  mean( abs(A2(:)) ) ]
  % [ std( A2(:) ) std( abs(A2(:)) ) ]
  
  set(la,'string','select a signal-void region outside the object');

  t2 = ones(length(t),1)*ginput(1)+[ t' t'];
  
  A2 = A( round(t2(:,2)), round(t2(:,1)) );
  ovrlay( round(t2(:,2)), round(t2(:,1)) ) = 0.7;
  set(ihdl,'cdata',ovrlay,'alphadata',0.7*ovrlay,'tag','overlay');
  
  cnr = vec( round([ t2(1,:) ; t2(size(t2,1),:) ]) + [ -0.5 -0.5; 0.5 0.5 ] );
  
  axis(get(ihdl,'parent'));
  l2 = line([ cnr(1) cnr(1) cnr(2) cnr(2) cnr(1) ], ... 
            [ cnr(3) cnr(4) cnr(4) cnr(3) cnr(3) ], ... 
            'color',[ 1 1 1 ], 'tag', 'box2');
  

  % [ mean( abs(A2(:)) )/1.253 std( abs(A2(:)) )/0.665 ]
  set(la,'string','');

  mhdl = findobj(fig,'tag','mhdl');
  curmde = get(mhdl,'value');
  shdl = findobj(fig,'tag','shdl');
  curval = get(shdl,'value');

  if isempty(curval),
    idx = 1;
  else,
    idx = get(shdl,'min'):get(shdl,'max');
  end;

  set(mhdl,'value',1); fprintf('\n');
  for cnt=1:length(idx),
    if ~isempty(curval),
      set(shdl,'value',cnt);
    end;
    eval([ mfilename '(''redraw'',' num2str(fig) ')' ]);
    A = get(ia,'cdata'); % field( get(fig,'userdata'), 'im_X' );
    A2 = A( round(t1(:,2)), round(t1(:,1)) );
    s = mean( abs(A2(:)) );
    
    A2 = A( round(t2(:,2)), round(t2(:,1)) );
    fprintf('Image %d - measured SNR: %f\n', cnt, s / std( abs(A2(:)) )/0.665 );
    pause(0.5);
  end;
  fprintf('\n');
  set(mhdl,'value',curmde);
  set(shdl,'value',curval);
  eval([ mfilename '(''redraw'',' num2str(fig) ')' ]);

elseif strcmp( in1, 'printall'),

  if ischar(parm1), 
    fname = parm1; 
  else
    fname = input('enter an base output filename:','s');
  end;
  
  sz = getfield( get(fig,'userdata'), 'size_X' );
  if length(sz)<3, sz(3) = 1; end;
  % if length(sz) > 2; ni = sz(3); else, ni = 1; end;
  shdl = getfield( get(fig,'userdata'), 'shdl');
  if length(shdl) == 2,
    for cnt=1:sz(3),
      set( shdl(1), 'value', cnt );
      for cnt2=1:sz(4),
        set( shdl(2), 'value', cnt2 );
        eval([ mfilename '(''redraw'',' num2str(fig) ')' ]);
        pause(0.1);
        
        cb = [ mfilename '( ''print'', ' num2str(fig) ',' ...
               sprintf(' ''%s_r%02d_i%02d.png'' ',fname,cnt,cnt2) ')' ];
        disp( cb );
        eval( cb );
      end;
    end;
  else
    for cnt=1:sz(3),
      if sz(3)~=1, set( shdl(length(shdl)), 'value', cnt ); end;

      eval([ mfilename '(''redraw'',' num2str(fig) ')' ]);
      pause(0.1);
        
      cb = [ mfilename '( ''print'', ' num2str(fig) ',' ...
             sprintf(' ''%s_%03d.png'' ',fname,cnt) ')' ];
      disp( cb );
      eval( cb );
    end;
  end;
  
elseif strcmp( in1, 'print'),

  if ischar(parm1), 
    fname = parm1; 
  else
    fname = input('enter an output filename:','s');
  end;
  figure(fig);

  mhdl = findobj(fig,'tag','mhdl');
  if isempty(mhdl),
    rsp = 1; 
  else,
    rsp = get( mhdl, 'value' );
  end;

  TTT = get(gca,'children');
  if length(TTT) ~= 1,
    TTT = findobj(TTT,'type','image');
  end;
  im_X = get(TTT,'CData');

  C = colormap;
  switch rsp,
   case 1,
    TTT = abs( im_X );  ttl = 'Intensity values';
    maxca = 'mxt'; minca = 'mnt';
   case 2,
    TTT = ( im_X ); ttl = 'Phase values';
    maxca = 'mxp'; minca = 'mnp';
   case 3,
    TTT = real( im_X );  ttl = 'Real values';
    maxca = 'mxr'; minca = 'mnr';
   case 4,
    TTT = imag( im_X );  ttl = 'Imaginary values';
    maxca = 'mxi'; minca = 'mni';
   case 5,
    TTT = ( im_X ); ttl = 'GOP';
  end;

  AAA = [];
  if (rsp~=5) 
    if iscell(TTT)
      AAA = get( findobj(fig,'tag','image mask'), 'AlphaData' );
      TTT = TTT{2};
    end;
    TTT = [ TTT-min(TTT(:)) ] * size(C,1) / max(TTT(:));
  
    cmx = max(TTT(:));
    as = round(axis);
    as([ 4 2 ]) = min([ as([4 2]); size(TTT) ]);
    as([ 1 3 ]) = max([ as([1 3]); [ 1 1 ] ]);

    TTT = TTT( as(3):as(4), as(1):as(2) );

    % TTT = TTT - min(TTT(:));
    if 1,
      if ( get( chdl, 'value' ) & ( rsp ~= 5 ) ),
        cmx = getfield( get(fig,'userdata'), maxca );
      elseif ( rsp == 5 ),
        cmx = size(C,1);
      else,
        cmx = max(TTT(:));
      end;
    end;
  
    TTT2 = floor(abs(TTT).*size(C,1)./cmx);
  else
    TTT2=TTT;
  end;

  if isempty(AAA),
    imwrite( TTT2, C, fname );
  else
    RGB = ind2rgb(TTT2,C);
    imwrite( RGB, fname, 'Alpha', 1-AAA );
  end;

elseif strcmp( in1, 'toggle'),

%   if isreal(im_X(parm1)),
%     im_X(parm1) = not( im_X(parm1) );
%   end;
% 
%   img_data = get(fig,'userdata');
%   if ~isempty(shdl),
%     img_data.im_X(:,:,round(get(shdl,'value'))) = im_X;
%   else
%     img_data.im_X = im_X;
%   end;
%   set( fig, 'userdata', img_data ); 
%   % setfield( set(1,'userdata'), 'im_X', img_data );
%   % getfield( get(1,'userdata'), 'im_X' )
%   im('redraw',fig);

elseif strcmp( in1, 'playmovie'),

  shdl = getfield( get(fig,'userdata'), 'shdl');
  mvhdl = getfield( get(fig,'userdata'), 'mvhdl');
  
  cb = get(mvhdl,'callback');
  set(mvhdl,'callback',[]);
  
  mx = get(shdl,'max');
  mn = get(shdl,'min');

  for cnt = mn:mx,
    if( get(mvhdl,'value') == 0 ), continue; end;

    set( shdl, 'value', cnt );
    im('redraw',fig);
    pause(0.5);
    M(cnt) = getframe(gca);
  end;
  
  while( get(mvhdl,'value') )
    movie(M,5,15);
  end;
  
  set(mvhdl,'callback',cb);
  set(mvhdl,'value',0);
  
elseif strcmp(in1,'setdisp'),

  mhdl = getfield( get(fig,'userdata'), 'menu');
  set(mhdl(:),'Checked','off');
  set(mhdl(parm1+1),'Checked','on');

  mhdl = getfield( get(fig,'userdata'), 'mhdl');
  set(mhdl,'value',parm1);

  im('redraw',fig);
end;

if nargout > 0, result = fig; end;

% Copyright (c) 2008 by W Scott Hoge (shoge -at- ieee -dot- org)
% Distributed under the terms of the NCIGT Fast Imaging Library 

function gopcoltable=goptab() 
% 
% Generates a GOP colortable 
% 
% by matsa & knutte 
% Dept. of Biomedical Engineering 
% Linköping University, Sweden 
% matsa@imt.liu.se  knutte@imt.liu.se 
% 
gopcoltable=[ 
93 234 63 
91 233 66 
89 232 69 
86 231 73 
84 230 76 
82 228 80 
80 227 83 
78 226 87 
76 225 91 
74 223 95 
72 222 99 
70 220 103 
69 219 107 
67 217 111 
65 216 115 
64 214 119 
62 212 123 
61 210 128 
59 208 132 
58 206 136 
57 204 141 
56 202 145 
55 200 149 
54 197 154 
54 195 158 
53 193 162 
53 190 166 
53 187 171 
53 185 175 
53 182 179 
53 179 183 
54 176 187 
55 173 191 
55 170 195 
56 167 199 
57 164 202 
59 161 206 
60 158 209 
61 155 212 
63 151 216 
64 148 219 
66 145 222 
68 142 225 
70 138 227 
72 135 230 
73 132 232 
75 128 235 
77 125 237 
79 121 239 
82 118 241 
84 115 243 
86 112 244 
88 108 246 
90 105 247 
92 102 249 
95 99 250 
97 95 251 
99 92 252 
101 89 253 
104 86 253 
106 83 254 
108 80 254 
111 77 255 
113 74 255 
115 71 255 
118 68 255 
120 65 255 
122 63 255 
125 60 254 
127 57 254 
130 54 253 
132 52 252 
135 49 252 
137 47 251 
140 44 250 
142 42 248 
145 39 247 
147 37 246 
150 34 244 
152 32 243 
155 30 241 
158 28 239 
160 26 237 
163 23 235 
165 21 233 
168 19 230 
171 17 228 
173 16 225 
176 14 222 
178 12 220 
181 10 217 
183 9 214 
186 7 210 
188 5 207 
191 4 204 
193 3 200 
195 1 196 
198 0 193 
200 0 189 
202 0 185 
204 0 181 
206 0 177 
208 0 173 
210 0 169 
212 0 164 
214 0 160 
215 0 156 
217 0 151 
218 0 147 
220 0 142 
221 0 138 
222 0 133 
223 0 129 
225 0 124 
226 0 120 
227 0 116 
227 0 111 
228 0 107 
229 0 103 
229 0 99 
230 0 94 
230 0 90 
231 0 86 
231 0 82 
231 0 79 
232 0 75 
232 0 71 
232 0 68 
232 0 64 
232 1 61 
232 2 58 
232 4 55 
232 5 51 
232 7 49 
232 8 46 
232 10 43 
232 11 40 
232 13 38 
232 15 36 
232 16 33 
232 18 31 
232 20 29 
232 22 27 
232 24 25 
232 26 23 
232 28 22 
232 30 20 
232 33 19 
232 35 17 
233 37 16 
233 39 14 
233 42 13 
233 44 12 
234 47 11 
234 50 10 
235 52 9 
235 55 8 
235 58 7 
236 61 7 
237 63 6 
237 66 5 
238 70 5 
238 73 4 
239 76 4 
240 79 3 
241 82 3 
241 86 2 
242 89 2 
243 93 2 
244 96 2 
245 100 1 
245 104 1 
246 108 1 
247 111 1 
248 115 1 
249 119 0 
250 123 0 
250 128 0 
251 132 0 
252 136 0 
252 140 0 
253 144 0 
254 149 0 
254 153 0 
254 157 0 
255 161 0 
255 166 0 
255 170 0 
255 174 0 
255 178 0 
255 183 0 
254 187 0 
254 191 0 
253 195 0 
253 199 0 
252 203 0 
251 206 0 
250 210 0 
249 213 0 
247 217 0 
246 220 0 
244 223 0 
242 226 0 
240 229 0 
238 232 0 
236 234 0 
234 237 0 
231 239 0 
229 241 0 
226 243 1 
224 245 1 
221 246 1 
218 248 1 
215 249 1 
212 250 1 
209 251 2 
206 252 2 
203 253 2 
200 253 3 
197 254 3 
193 254 4 
190 255 4 
187 255 5 
184 255 5 
181 255 6 
177 255 7 
174 255 7 
171 255 8 
168 254 9 
165 254 10 
162 254 11 
159 253 12 
156 253 13 
153 252 14 
150 252 15 
147 251 17 
144 250 18 
141 250 20 
138 249 21 
135 248 23 
132 248 24 
130 247 26 
127 246 28 
124 245 30 
122 245 32 
119 244 34 
117 243 36 
114 242 39 
112 241 41 
109 240 43 
107 240 46 
104 239 48 
102 238 51 
100 237 54 
97 236 57 
95 235 60 
]/255; 
