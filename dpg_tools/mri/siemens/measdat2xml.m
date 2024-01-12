function xml = measdat2xml(arg1);

mode = 0;
if isstruct(arg1)
  mode = 1;
  N = length(arg1);
elseif iscell(arg1)
  mode = 2;
  N = length(arg1);
else
  N = 1;
end;

for cntF=1:N,

  if (mode == 1)
    fname = arg1(cntF).name;
  elseif (mode == 2 )
    fname = arg1{cntF};
  else
    fname = arg1(cntF,:);
  end;
  fprintf('%s\n',fname);

fid = fopen(fname,'r');

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

  if (NScans>1)
    % loop through the scans, to select the data header for the last scan:
    cPos = datStart; 
    nscan_flag = NScans-1;
    while (nscan_flag); 
      fseek(fid,cPos,'bof'); 
      lngth = fread(fid, 1,'ubit25'); % fprintf('%d ',lngth); 
      cPos=cPos+lngth; 
      if (lngth< 1280)                     % reached the end of a scan
        % jump to next full 512 bytes
        if (lngth==0) cPos = cPos + 512 - mod(cPos,512); end;

        fseek(fid,cPos,'bof');
        hdrLength  = fread(fid,1,'uint32');
        nscan_flag = nscan_flag-1;
        if (nscan_flag)                 % if not the last scan, then
          cPos = cPos + hdrLength;      % skip the header, and move to the next data block
        end;
      end; 
    end;
  end;


else
  % VB style files
  data_version = 'VB';

  hdrLength = tmp-4;
  datStart = tmp;
end;

bytes_to_skip_hdr = datStart;

asc = fread(fid,bytes_to_skip_hdr,'char');
asc = char(asc');

fclose(fid);


for cnt=1:20; 
  if ( regexp( asc(cnt+(0:15)), '<XProtocol>' ) )
    break;
  end;
end;
asc = asc(cnt:end);

ndx = regexp( asc( (-40:0)+end ), '}' );
if ~isempty(ndx)
    asc = asc( 1: (end-40+max(ndx)) );
end;


fname2 = regexprep( fname,'.dat','.xml');


if ( strcmp(fname,fname2) )
  
else
  % keyboard;
  if nargout == 0
    fid = fopen(fname2,'w');
    fwrite(fid,asc,'char');
    fclose(fid);
  else
    xml = sprintf('%c',asc);
  end;
end;

end;

