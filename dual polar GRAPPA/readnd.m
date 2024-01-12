function [A,meta] = readnd(fname);
%
% READND - read a binary N-dimensional data file
%
% The function reads N-dimensional data saved by either the C or Matlab
% form of savend.
%
%    See Also SAVEND
%

%
% Copyright 2004-2008 W Scott Hoge (shoge at ieee dot org) 
% Distributed under the terms of the NCIGT Fast Imaging Library 


fid = fopen(fname,'rb');
if fid <0, error(['Couldn''t open ' fname]); end;

fseek(fid,0,'eof');
maxfln = ftell(fid);
fseek(fid,0,'bof');

fln = fread(fid,1,'integer*4');
hdr = fscanf(fid,'%c',4);          
if ~strcmp(hdr,'nddf')
  warning(['This does not appear to be a .nd file.  The file type field is ' ...
           'not ''nddf'' ']);
  return;
end;
 
szln = fread(fid,1,'integer*4');    
  sz = fread(fid,szln,'integer*4');
 hdr = fscanf(fid,'%c',4);

if ( strcmp(hdr(4),'c') ),
  cmplx = 2;
else
  cmplx = 1;
end;
   
if strcmp( hdr(1:3), 'dbl' ),
  data = fread(fid,cmplx*prod(sz),'double');
elseif strcmp( hdr(1:3), 'flt' ),
  data = fread(fid,cmplx*prod(sz),'float');
elseif strcmp( hdr(1:3), 'int' ),
  data = fread(fid,cmplx*prod(sz),'int16');
elseif strcmp( hdr(1:3), 'cha' ),
  data = fread(fid,cmplx*prod(sz),'uint8');
end;
fprintf('format: %s\n',hdr);

if ( strcmp(hdr(4),'c') ),
  data = data(1:2:length(data)) + j * data(2:2:length(data));
end;
if ( prod(size(data)) ~= prod(sz) ), 
  fprintf('declared: %10d\n',prod(sz));
  fprintf('    read: %10d\n',length(data));
  warning('Error: size of read data doesn''t match declared size');
else
  if length(sz)==1, sz = size(data)'; end;
  A = reshape( data, sz' );
end;

if ( ftell(fid) ~= maxfln ), 

  metaln = fread(fid,1,'integer*4');
  hdr = fscanf(fid,'%c',4);  
  
  if strcmp( hdr, 'meta' ),
    meta = fread(fid,metaln,'char');
  else,
    warning('Warning: number of bytes read from file is suspicious'); 
    disp([ ftell(fid) fln ]);
    keyboard;
  end;
end;

fclose(fid); 

