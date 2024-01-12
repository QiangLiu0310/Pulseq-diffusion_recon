function refdatgui(in1,parm1,parm2)

% refdatgui( k, [ lp, cp ] )

if ~ischar(in1)

  k  = in1;
  lp = []; cp = [];
  if nargin > 1, lp = parm1; end;
  if nargin > 2, cp = parm2; end;

  fig = figure;
  if ~isreal(fig) % check if 'fig' is a matlab.ui.Figure class object
     fig = fig.Number;
  end;

  cb = [ mfilename '(''redraw'', ', num2str(fig), ')' ];

  shdl = uicontrol('style','slider','position', [ 30 10 180 20 ], ...
                   'max', size(k,1)/2, 'min', -size(k,1)/2, 'value', 0, ...
                   'tag', 'shdl', ...
                   'tooltip','click here to adjust the linear phase correction', ...
                   'sliderstep', [ 0.001 0.01 ], ...
                   'callback', cb );
  ehdl = uicontrol('style','edit','position', [ 210 10 50 20 ], ...
                   'value', 0, 'string', '0', ...
                   'tag', 'hdl1', 'callback', cb );

  cs = uicontrol('style','slider','position', [ 30 30 180 20 ], ...
                   'max', pi, 'min', -pi, 'value', 0, ...
                   'tag', 'shdl', ...
                   'tooltip','click here to adjust the constant phase correction', ...
                   'sliderstep', [ 0.001 0.01 ], ...
                   'callback', cb );
  ce = uicontrol('style','edit','position', [ 210 30 50 20 ], ...
                   'value', 0, 'string', '0', ...
                   'tag', 'hdl1', 'callback', cb );

  set(fig,'userdata',struct('k', double(in1), ...
                            'lp', lp, ...
                            'cp', cp, ...
                            'shdl',shdl, ...
                            'ehdl',ehdl, ...
                            'cs',cs, ...
                            'ce',ce  ...
                            ) );

  eval( cb );
  
elseif strcmp(in1,'redraw')

  if nargin==1,
    fig = get(gcbo,'parent');
  else
    fig = parm1;
  end;

  usrdat = get(fig,'userdata');
  
  %% linear phase GUI update
  if ( gcbo == usrdat.ehdl ),
    val = str2num( get(gcbo,'string') );
    set( usrdat.shdl, 'value', val );
  end;
  
  lphz = get( usrdat.shdl, 'value' );
  set( usrdat.ehdl, 'string', num2str(lphz) );

  %% constant phase GUI update
  if ( gcbo == usrdat.ce ),
    val = str2num( get(gcbo,'string') );
    set( usrdat.cs, 'value', val );
  end;
  
  cphz = get( usrdat.cs, 'value' );
  set( usrdat.ce, 'string', num2str(cphz) );

  k = usrdat.k;
  indx = [ 0 : (size(k,2)-1) ] - size(k,2)/2;
  % k = apply_readout_shift( k, [ cphz; lphz/size(k,1) ] );
  k = phzapply( k.', [ cphz; lphz ] ).';
  
  subplot(121);
  imagesc( abs(k) ); axis('image');
  
  subplot(122);
  imagesc( ffts(abs(ifft2(ffts(k,1:2))),1:2) ); axis('image');

end;

