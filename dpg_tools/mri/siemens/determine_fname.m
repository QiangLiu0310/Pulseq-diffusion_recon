MID = regexp(pwd,'MID\d+','match');
f = dir([ '../*' MID{1}  '*.dat']);
if isempty(f)
  error('  :-(  no meas.dat files found   :-(');
end;

fname = [ '../' f(1).name ];
clear f MID;

disp([ 'fname = ''' fname '''' ]);