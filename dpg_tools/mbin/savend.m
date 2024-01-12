function savend(fname,A,frmt,meta);
%
% savend(fname,A,[frmt]);
%
% SAVEND - save an N-dimensional variable using single precision
%
% The function wraps the output data with a simple header to facilitate
% automatic reading of N-dimensional data in both C and Matlab.
%
%   fname - the name of the output file
%       A - the variable to save
%    frmt - [optional] the precision.  single:'flt' (default) or double:'dbl'
%
%    See Also READND
%

% Copyright 2003-2004, W Scott Hoge (shoge at ieee dot org) 
% All rights reserved. Licensed according to the GNU General Public Licence
% (see http://www.gnu.org)
% 

if nargin<3, frmt='flt'; end;
if nargin<4, meta=''; end;

fid = fopen(fname,'wb');
if fid <0, error(['Couldn''t open ' fname]); end;

sz = size(A);

if strcmp(frmt,'int')
  sizeof = 2;
elseif strcmp(frmt,'flt')
  sizeof = 4;
else,
  sizeof = 8;
end;

total = length(A(:))*sizeof + 4*4 + 4*length(sz);
if ~isreal(A),
  total = total + length(A(:))*sizeof;
end;

fwrite(fid,total,'integer*4');
fprintf(fid,'nddf');
fwrite(fid,length(sz),'integer*4');
fwrite(fid,sz(:),'integer*4');

if strcmp(frmt,'int'),
  if isreal(A),
    fwrite(fid,'intr');
    fwrite(fid,A(:),'int16');
  else
    fwrite(fid,'intc');
    fwrite(fid,[ real(A(:))'; imag(A(:))' ],'int16');
  end;
elseif strcmp(frmt,'char'),
  if isreal(A),
    fwrite(fid,'char');
    fwrite(fid,A(:),'int8');
  else
    fwrite(fid,'char');
    fwrite(fid,abs(A(:)),'int8');
  end;
elseif strcmp(frmt,'flt'),
  if isreal(A),
    fwrite(fid,'fltr');
    fwrite(fid,A(:),'float');
  else
    fwrite(fid,'fltc');
    fwrite(fid,[ real(A(:))'; imag(A(:))' ],'float');
  end;
else,
  if isreal(A),
    fwrite(fid,'dblr');
    fwrite(fid,A(:),'double');
  else
    fwrite(fid,'dblc');
    fwrite(fid,[ real(A(:))'; imag(A(:))' ],'double');
  end;
end;

if ( length(meta)>0 )
  fwrite(fid,length(meta)+4,'int32');
  fwrite(fid,'meta');
  fwrite(fid,meta,'char');
end;

fclose(fid); 