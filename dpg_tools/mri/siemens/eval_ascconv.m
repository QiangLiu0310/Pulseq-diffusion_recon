function ASCCONV = eval_ascconv(inp)
%
% this function converts variables expressed in C-code style to a form
% that Matlab can interpret and then store them as a structured variable.
%

fid = fopen(inp,'r');
tmp = fread(fid,1,'uint32');

if ( tmp < 10000 )
  % VD style files
  data_version = 'VD';

  NScans = fread(fid,1,'uint32');
  measID = fread(fid,1,'uint32');
  fileID = fread(fid,1,'uint32');
  
  % measOffset: points to beginning of header, usually at 10240 bytes
  measOffset = fread(fid,1,'uint64');
  measLength = fread(fid,1,'uint64');
  fseek(fid,measOffset,'bof');
  hdrLength  = fread(fid,1,'uint32');
  datStart   = measOffset + hdrLength;
else
  % VB style files
  data_version = 'VB';

  hdrLength = tmp-4;
  datStart = tmp;
end;

bytes_to_skip_hdr = datStart;


asc = fread(fid,hdrLength,'char');
fclose(fid);
asc = char(asc');


t1 = strfind(asc,'### ASCCONV BEGIN ');
t2 = strfind(asc,'### ASCCONV END ###');

% extract the ASCCONV section, append each line w/ a semi-colon, and 
% add in a top variable name (here, 'ASCONV') that will be returned.
if (data_version == 'VD')
  tmp = regexprep( asc(t1(end):t2(end)-2), '\n', ';\nASCCONV.' ); 
else,
  tmp = regexprep( asc(t1(1):t2(1)-2), '\n', ';\nASCCONV.' ); 
end;

tmp = [ tmp sprintf(';\n') ];

tmp = regexprep( tmp,'\[(\d+)\]', '[$1+1]' ); % add '1' to each specified array index,
tmp = regexprep( tmp,'[','(');                % convert C-style array braces
tmp = regexprep( tmp,']',')');                % ... to Matlab-style
tmp = regexprep( tmp,'#','%');                % convert comment lines
tmp = regexprep( tmp,'(ASC\S*__attrib)','% $1');                % comment out  __attribute__ lines
tmp = regexprep( tmp,'"+','''');               % convert C-strings to Matlab-strings
tmp = regexprep( tmp,'\s0x([a-f\d]+)','hex2dec(''$1'');'); % convert all 0xDD numbers to hex2dec
tmp = regexprep( tmp,'\s0x([a-f\d]+) (% ''[A-Z]'');','hex2dec(''$1''); $2'); % convert all 0xDD numbers to hex2dec

indx = [1 regexp(tmp,'\n')];

eval(tmp);                        
% eval the string as a list of variable declarations

% for cnt=1:length( ASCCONV.sSliceArray.asSlice ),
%   if ~isfield( ASCCONV.sSliceArray.asSlice(cnt).sPosition, 'dSag'),
%     ASCCONV.sSliceArray.asSlice(cnt).sPosition.dSag = 0;
%   end;
%   if ~isfield( ASCCONV.sSliceArray.asSlice(cnt).sPosition, 'dCor'),
%     ASCCONV.sSliceArray.asSlice(cnt).sPosition.dCor = 0;
%   end;
%   if ~isfield( ASCCONV.sSliceArray.asSlice(cnt).sPosition, 'dTra'),
%     ASCCONV.sSliceArray.asSlice(cnt).sPosition.dTra = 0;
%   end;
% end;
