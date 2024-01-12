function hideobjs(in1,fig),
% HIDEOBJS prevent certain figure objects from printing by toggling the visible flag
%
  % usage:  hideobjs(0)       sets objects visibility to 'off'
  %         hideobjs(1)       sets objects visibility to 'on'
%
  types  = { 'popupmenu', 'radiobutton', 'checkbox', 'slider', 'pushbutton', 'text', 'edit' };
action = { 'off', 'on' };

if nargin == 0, in1 = 0; end;
if nargin < 2, fig = gcf; end;


for cnt1=1:length(types),
	   hdl = findobj( fig, 'style', types{cnt1} );
for cnt2=1:length(hdl),
	   set( hdl(cnt2), 'visible', action{ mod(in1,length(action))+1 } );
end;
end;
