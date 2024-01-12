function planevu( action, argv );


if ~ischar(action),

  if (nargin<2),
    fig = figure;
    if ~isreal(fig),
        fig = fig.Number;
    end;
  else,
    fig = argv;
    if ~isreal(fig),
        fig = fig.Number;
    end;
    figure(fig);
    clf;
  end;

  in = action;
  if isreal(in)

  else
%    in = abs(in);
  end;


  sz = size(in);

  pos = round([ sz(1:3)/2 ]);

  Im = squeeze( abs(in(:,:,:,1)) );

  ahdl(1) = subplot(221);
  ihdl(1) = image( (Im( :, :, round(end/2), 1 )) );
  xhdl(1) = xlabel(['-  ' num2str(pos(1)) '  +']);
  set(xhdl(1),'userdata', [ pos(1) ]);

  ahdl(2) = subplot(222);
  ihdl(1) = image( squeeze((Im( :, round(end/2), :, 1 ))) );
  xhdl(2) = xlabel(['-  ' num2str(pos(2)) '  +']);
  set(xhdl(2),'userdata', [ pos(2) ]);

  ahdl(3) = subplot(223);
  ihdl(3) = image( squeeze((Im( round(end/2), :, :, 1 ))) );
  xhdl(3) = xlabel(['-  ' num2str(pos(3)) '  +']);
  set(xhdl(3),'userdata', [ pos(3) ]);

  set(ahdl(1),'tag','axis1');
  set(ahdl(2),'tag','axis2');
  set(ahdl(3),'tag','axis3');

  ahdl(4) = subplot(224);
  if (length(sz)>3),
    set(ahdl(4),'visible','off');
  end;
  set(fig,'units','pixels');
  tmp = get(fig,'pos');

  axscl=(tmp(4)-44-60-60)/(sz(1)+sz(3)); % 4;  5.16;
  set(ahdl(3),'units','pixels','pos',[ 100 60 axscl*[sz(2) sz(1)] ]);
  set(ahdl(1),'units','pixels','pos',[ 100 60+60+axscl*sz(1) axscl*[sz(2) sz(3)] ]);
  set(ahdl(2),'units','pixels','pos',[ 100+axscl*sz(2)+60 60+60+axscl*sz(1) axscl*[sz(1) sz(3)] ]);

  cb = [ mfilename '(''redraw'', ', num2str(fig), ')' ];

  if (length(sz)>3),
    % ahdl(4) = asubplot(2,2,4,'units','pixels','pos',[ tmp(3) 60 axscl*[sz(2) sz(3)] ]);

    
    phdl = plot( (squeeze( Im(pos(1),pos(2),pos(3),:,1) )) );
    set(phdl,'tag','tplot','hittest','off'); hold on;

    phdl = plot( (squeeze( Im(pos(1),pos(2)+floor(size(Im,2)/2),pos(3),:,1) )),'k');
    set(phdl,'tag','tplot2','hittest','off','visible','off'); hold on;

    phdl = plot( 1, (squeeze( Im(pos(1),pos(2),pos(3),1) )), 'ok', ...
                 'markersize', 8 );
    set(phdl,'tag','tpnt','hittest','off');
    hold off;
    tmp = axis;
    axis([ 1 sz(4) tmp(3:4) ]);
    xhdl(4) = xlabel(['-  ' num2str(1) '  +']);
    set(xhdl(4),'userdata', 1 );

    set(ahdl(4),'tag','axis4');
    set(ahdl(4),'units','pixels','pos',[ 100+axscl*sz(2)+60 60 axscl*[sz(1) sz(1)] ]);

    p2hdl = uicontrol( 'style', 'checkbox', 'string', 'plot2', ...
                       'value', 0, 'tag', 'plot2', ...
                       'position', [ 90 00 80 20 ], ...
                       'tooltip', ['plots 2 pixels in the fourth sub-plot'], ...
                       'callback', cb );

  end;
  set(ahdl,'buttondownfcn',@onclick);
  set(xhdl,'buttondownfcn',@increment );

  uicontrol( 'style', 'checkbox', 'string', 'autoscale', ...
             'value', 1, 'tag', 'imgscl', ...
             'position', [ 10 20 80 20 ], ...
             'tooltip', ['switches image display between ''image'' and ' ...
                      '''imagesc''' ], ...
             'callback', cb );

  uicontrol( 'style', 'checkbox', 'string', 'cmpr a:b', ...
             'value', 0, 'tag', 'cmprab', ...
             'position', [ 10 40 80 20 ], ...
             'tooltip', ['display image as two 4D side-by-side volumes'], ...
             'callback', cb );

  uicontrol( 'style', 'pushbutton', 'string', 'output', ...
             'tag', 'print', ...
             'position', [ 10 60 80 20 ], ...
             'tooltip', ['output .png versions of each view'], ...
             'callback',  [ mfilename '(''print'', ', num2str(fig), ')' ] );

  mhdl = uicontrol( 'style', 'checkbox', 'string', 'mag/phz', ...
                    'value', 0, 'tag', 'imgphz', ...
                    'position', [ 10 0 80 20 ], ...
                    'tooltip', ['switches between magnitude and phase values of pixels '], ...
                    'callback', cb );
  if (isreal(in))
    set(mhdl,'visible','off');
  end;


%   uicontrol( 'style', 'checkbox', 'string', 'autoaxis', ...
%              'value', 1, 'tag', 'imgaxs', ...
%              'position', [ 10 0 80 20 ], ...
%              'tooltip', ['switches axis ''auto'' and ' ...
%                       '''normal''' ], ...
%              'callback', cb );

  udat = struct( 'sz', sz, ...
                 'im', in, ...
                 'mh', mhdl, ...
                 'ah', ahdl );
  set(fig,'userdata',udat);

  action = 'redraw';

  set(fig, 'Name', [ mfilename '( ' inputname(1) ' )' ]);

  clear('in');
else,
  fig = argv;
end;

if strcmp(action,'redraw'),

  ahdl = getfield( get(fig,'userdata'), 'ah' );

%   str = get(get(ahdl(1),'xlabel'), 'string' );
%   pos(1) = str2num( str(regexp(str,'\w')) );
%   str = get(get(ahdl(2),'xlabel'), 'string' );
%   pos(2) = str2num( str(regexp(str,'\w')) );
%   str = get(get(ahdl(3),'xlabel'), 'string' );
%   pos(3) = str2num( str(regexp(str,'\w')) );
  pos(1) = get(get(ahdl(1),'xlabel'), 'userdata' );
  pos(2) = get(get(ahdl(2),'xlabel'), 'userdata' );
  pos(3) = get(get(ahdl(3),'xlabel'), 'userdata' );

  sz = getfield( get(fig,'userdata'), 'sz' );

  if ( length(sz)>3),
    pt = get(get(ahdl(4),'xlabel'), 'userdata' );
  else,
    pt = 1;
  end
  if isempty(pt)
    pt = 1;
  end;

  if ( get( findobj(fig,'tag','imgphz'), 'value' ) == 1);
    Im = angle(getfield( get(fig,'userdata'), 'im' ));
  else
    Im = abs(getfield( get(fig,'userdata'), 'im' ));
  end;
  if ~isfloat(Im)
    Im = double(Im);
  end;
  
  pos( pos < 1 ) = 1;
  indx = pos > sz(1:3);
  if length(indx)>0,
    pos( indx ) = sz( indx );
  else,

  end;

  scl  = get( findobj(fig,'tag','imgscl'), 'value' );
  cmpr = get( findobj(fig,'tag','cmprab'), 'value' );

  maxval = max(abs(Im(:)));
  minval = min(abs(Im(:)));
  
  sclfac = length(colormap) / ( max((vec(Im(:,:,:,pt)))) - min((vec(Im(:,:,:,pt)))) );
  mnfac  = min((vec(Im(:,:,:,pt))));

  sclfac1 = length(colormap) / ( max((vec(Im(pos(1),:,:,pt)))) - min((vec(Im(pos(1),:,:,pt)))) );
  sclfac2 = length(colormap) / ( max((vec(Im(:,pos(2),:,pt)))) - min((vec(Im(:,pos(2),:,pt)))) );
  sclfac3 = length(colormap) / ( max((vec(Im(:,:,pos(3),pt)))) - min((vec(Im(:,:,pos(3),pt)))) );
  mnfac1 = min((vec(Im(pos(1),:,:,pt))));
  mnfac2 = min((vec(Im(:,pos(2),:,pt))));
  mnfac3 = min((vec(Im(:,:,pos(3),pt))));

  autoaxis = 1; % get( findobj(fig,'tag','imgaxs'), 'value' );

  figure(fig);

  tmp = get(fig,'pos');
  axscl=(tmp(4)-44-60-60)/(sz(1)+sz(3));
  set(ahdl(3),'units','pixels','pos',[ 100 60 axscl*[sz(2) sz(1)] ]);
  set(ahdl(1),'units','pixels','pos',[ 100 60+60+axscl*sz(1) axscl*[sz(2) sz(3)] ]);
  set(ahdl(2),'units','pixels','pos',[ 100+axscl*sz(2)+60 60+60+axscl*sz(1) axscl*[sz(1) sz(3)] ]);
  if (length(ahdl)>3),
    set(ahdl(4),'units','pixels','pos',[ 100+axscl*sz(2)+60 60 axscl*[sz(1) sz(1)] ]);
  end;
  set(fig,'pos',[tmp(1:2) 100+60+60+axscl*[sz(2)+sz(1)] tmp(4)]);

  subplot(ahdl(1));
  ihdl(1) = findobj( get(ahdl(1),'children'), 'type', 'image' );
  delete( findobj( get(ahdl(1),'children'), 'type', 'line' ) );

  tmpdat = (squeeze(Im(pos(1),:,:,pt))).';
  if scl,
    set( ihdl(1), 'CData', ( tmpdat - mnfac1) * sclfac1 );
  else,
    set( ihdl(1), 'CData', ( tmpdat - mnfac) * sclfac );
  end;
  if autoaxis,
    axis('image');
  else
    axis('normal');
  end;
  set( get(ahdl(1),'xlabel'), ...
       'string', [ '-  ' num2str( pos(1) ) '  +' ], ...
       'userdata', pos(1) ...
       );

  l1 = line([ pos(2) pos(2); 1 sz(2) ]', [ 1 sz(3); pos(3) pos(3) ]' );
  set(l1,'color',[1 1 1]);

  % the 'axial' view:
  subplot(ahdl(2));
  ihdl(2) = findobj( get(ahdl(2),'children'), 'type', 'image' );
  delete( findobj( get(ahdl(2),'children'), 'type', 'line' ) );

  if cmpr
      pos2 = pos(2) + round(sz(2)/2);
      if (pos2 > sz(2)), pos2 = pos2 - sz(2); end;
      if (pos2 < 1), pos2 = pos2 + sz(2); end;
      tmpdat = [ (squeeze(Im(:,pos(2),:,pt))).' (squeeze(Im(:,pos2,:,pt))).'  ];
  else
      tmpdat = (squeeze(Im(:,pos(2),:,pt))).';
  end;
  
  if scl,
    set( ihdl(2), 'CData', ( tmpdat - mnfac2 ) * sclfac2 );
  else,
    set( ihdl(2), 'CData', ( tmpdat - mnfac ) * sclfac );
  end;
  if autoaxis,
    axis('image');
  else,
    axis('normal');
  end;
  set( get(ahdl(2),'xlabel'), ...
       'string', [ '-  ' num2str( pos(2) ) '  +' ], ...
       'userdata', pos(2) ...
       );
  l2 = line([ pos(1) pos(1); 1 sz(1) ]', [ 1 sz(3); pos(3) pos(3) ]' );
  set(l2,'color',[1 1 1]);

  % the 'axial' view
  subplot(ahdl(3));
  ihdl(3) = findobj( get(ahdl(3),'children'), 'type', 'image' );
  delete( findobj( get(ahdl(3),'children'), 'type', 'line' ) );
  if scl,
    set( ihdl(3), 'CData', ((squeeze(Im(:,:,pos(3),pt))) - mnfac3 ) * sclfac3 );
  else,
    set( ihdl(3), 'CData', ((squeeze(Im(:,:,pos(3),pt))) - mnfac ) * sclfac );
  end;
  if autoaxis,
    axis('image');
  else
    axis('normal');
  end;
  set( get(ahdl(3),'xlabel'), ...
       'string', [ '-  ' num2str( pos(3) ) '  +' ], ...
       'userdata', pos(3) ...
       );
  l3 = line(  [ 1 sz(2); pos(2) pos(2) ]', [ pos(1) pos(1); 1 sz(1) ]' );
  set(l3,'color',[1 1 1]);

  if (length(sz)>3),
    
    if ( get( findobj(fig,'tag','plot2'), 'value' ) == 1 ) 
      subplot(ahdl(3));
      l3 = line(  mod([ pos(2)+sz(2)/2 pos(2)+sz(2)/2 ]'-1,sz(2))+1, [ 1 sz(1) ]' );
      set(l3,'color',[1 1 1]);
    end;

    a4 = ahdl(4);
    
    if 1,
      v = vec(squeeze( Im(pos(1),pos(2),pos(3),:) ));
    else,
      v = vec(squeeze( Im(pos(1),pos(2),pos(3),:,:) ));
    end;

    phdl = findobj( get(a4,'children'), 'tag', 'tplot' );
    set(phdl, 'xdata', 1:length(v), 'ydata', v );

    if ( get( findobj(fig,'tag','plot2'), 'value' ) == 1 ) 
      if 1,
        v2 = vec(squeeze( Im(pos(1), mod( (pos(2)+size(Im,2)/2)-1,size(Im,2))+1,pos(3),:) ));
      else,
        v2 = vec(squeeze( Im(pos(1),pos(2),pos(3),:,:) ));
      end;

      phdl = findobj( get(a4,'children'), 'tag', 'tplot2' );
      set(phdl, 'ydata', v2, 'visible', 'on' );
    else,
      phdl = findobj( get(a4,'children'), 'tag', 'tplot2' );
      set(phdl, 'visible', 'off' );
    end;
    
    phdl = findobj( get(a4,'children'), 'tag', 'tpnt' );
    set(phdl, 'xdata', pt, 'ydata', v(pt) );
    v(pt);

    axis(ahdl(4),'auto');
    tmp = axis(ahdl(4));

    if scl,
      axis(ahdl(4),[ 1 sz(4) tmp(3:4) ]);
    else,
      if ( get( findobj(fig,'tag','imgphz'), 'value' ) == 1);
        axis(ahdl(4),[ 1 sz(4) -pi pi ]);
      else
        axis(ahdl(4),[ 1 sz(4) minval maxval ]);
      end;
    end;

    y = get( get(ahdl(4),'xlabel'), 'userdata');  
    set( get(ahdl(4),'xlabel'), 'string',['-   ' num2str(y) '   +'],'userdata',y);  
  else
    set(ahdl(4),'visible','off');
  end;

  set(ihdl,'hittest','off');
  set(ahdl(1:3),'buttondownfcn',@onclick);
  set(ahdl(1),'tag','axis1');
  set(ahdl(2),'tag','axis2');
  set(ahdl(3),'tag','axis3');

elseif strcmp(action,'print'),

    ahdl = getfield( get(fig,'userdata'), 'ah' );
    cmap = colormap;
    lbls = {'sa','co','ax'};
    for cnt=1:3
      ch = get(ahdl(cnt),'children');
      img = findobj( ch, 'type','image');
      I = get(img,'CData');
      pt = get(get(ahdl(cnt),'xlabel'),'userdata');
      imwrite( I, cmap, sprintf('img_3pv_%s%03d.png',lbls{cnt},pt) );
    end;

end;

return;

function increment(gcbo,eventdata,handles)

% keyboard;
ca = get(gcbo,'parent');
curpt = get(ca,'currentpoint');
fig = get(ca,'parent');

sz = getfield( get(fig,'userdata'), 'sz' );

axs = axis(ca);

switch( get(ca,'tag') )
 case 'axis1'
  maxpos = sz(1);
 case 'axis2'
  maxpos = sz(2);
 case 'axis3'
  maxpos = sz(3);
 case 'axis4'
  maxpos = sz(4);
end;
% str = get(gcbo,'string');
% curpos = str2num( star( regexp(str,'\w') ) );
curpos = get(gcbo,'userdata');

if ( curpt(1) > (axs(2)-axs(1))/2 + axs(1) ),
  if (curpos < maxpos),
    curpos = curpos+1;
  end;
else,
  if (curpos > 1),
    curpos = curpos-1;
  end;
end;

set(gcbo,'string',['-   ' num2str(curpos) '   +'],'userdata',curpos);


planevu('redraw',fig);

return;

function onclick(gcbo,eventdata,handles)

% http://blogs.mathworks.com/pick/2007/12/26/advanced-matlab-buttondownfcn/
curpt = get(gcbo,'currentpoint');

fig = get(gcbo,'parent');

x = round(curpt(3));
y = round(curpt(1));
if x<1; x=1; end;
if y<1; y=1; end;
% fprintf('%d %d\n',x,y);
ahdl = getfield( get(fig, 'userdata'), 'ah' );

switch( get(gca,'tag') )
 case 'axis1'
  pos(2) = y; pos(3) = x;
  set( get(ahdl(3),'xlabel'), 'string', num2str( pos(3) ),'userdata',pos(3));
  set( get(ahdl(2),'xlabel'), 'string', num2str( pos(2) ),'userdata',pos(2) );
 case 'axis2'
  pos(3) = x; pos(1) = y;
  set( get(ahdl(1),'xlabel'), 'string', num2str( pos(1) ),'userdata',pos(1) );
  set( get(ahdl(3),'xlabel'), 'string', num2str( pos(3) ),'userdata',pos(3) );
 case 'axis3'
  pos(1) = x; pos(2) = y;
  set( get(ahdl(2),'xlabel'), 'string', num2str( pos(2) ),'userdata',pos(2) );
  set( get(ahdl(1),'xlabel'), 'string', num2str( pos(1) ),'userdata',pos(1) );
 case 'axis4'
  pos(1) = x; pos(2) = y;
  set( get(ahdl(4),'xlabel'), 'string',['-   ' num2str(y) '   +'],'userdata',y);

end;

planevu('redraw',fig);
