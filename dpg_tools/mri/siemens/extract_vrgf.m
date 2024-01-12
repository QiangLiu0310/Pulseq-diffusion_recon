function vrgf = extract_vrgf( filename )

fid = fopen(filename,'r');
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

  cnt = 1;
  if (NScans>1)
    % loop through the scans, to select the data header for the last scan:
    cPos = datStart; 
    nscan_flag = NScans-1;
    while (nscan_flag); 
      fseek(fid,cPos,'bof'); 
      lngth = fread(fid, 1,'ubit25'); 
      cPos=cPos+lngth;
      % fprintf('%4d cPos: %8d  lngth: %8d\n',cnt,cPos,lngth);
      cnt=cnt+1;
      if (lngth == 0 ) % < 1280)    %                 % reached the end of a scan
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
  NScans = 1;
  hdrLength = tmp;
  datStart = tmp;
end;

% bytes_to_skip_hdr = datStart;
asc = fread(fid,hdrLength-4,'char');
asc = char(asc');
% fid = fopen('acs.tmp','w'); fwrite(fid,asc,'char'); fclose(fid);
vrgf = regrid_epi_rampsamp( asc );

return;

function vrgf = regrid_epi_rampsamp( asc )

vrgf  = 1; 
seflg = 0;

tmp  = '<ParamLong."FlattopTime">  {';
indx = strfind(asc,tmp);
offset = 0;
if isempty(indx), 
  tmp  = '<ParamLong."alRegridFlattopTime">';
  indx = strfind(asc,tmp);
  % keyboard;
  if isempty(indx), 
    return; 
  else
    seflg = 1;
    offset = 10;
  end;
end;
t_ft = sscanf( asc( indx(1)+length(tmp)+(0:32)+offset ), '%d' );
t_ft = t_ft(1);

if (seflg)
  tmp  = '<ParamLong."alRegridRampupTime">';
else
  tmp  = '<ParamLong."RampupTime">  {';
end;
indx = strfind(asc,tmp);
t_ru = sscanf( asc( indx(1)+length(tmp)+(0:20)+offset ), '%d' );
t_ru = t_ru(1);

if (seflg)
  tmp  = '<ParamLong."alRegridRampdownTime">';
else
  tmp  = '<ParamLong."RampdownTime">  {';
end;
indx = strfind(asc,tmp);
t_rd = sscanf( asc( indx(1)+length(tmp)+(0:20)+offset ), '%d' );
t_rd = t_rd(1);

if (seflg)
  tmp  = '<ParamLong."alRegridDelaySamplesTime">';
else
  tmp  = '<ParamLong."DelaySamplesTime">  {';
end;
indx = strfind(asc,tmp);
t_ds = sscanf( asc( indx(1)+length(tmp)+(0:20)+offset ), '%d' );
t_ds = t_ds(1);

if (seflg)
  tmp  = '<ParamLong."alRegridDestSamples">';
  offset = 10;
else
  tmp  = '<ParamLong."DestSamples">  {';
  offset = 0;
end;
indx = strfind(asc,tmp);
inres = sscanf( asc( indx(1)+length(tmp)+(0:20)+offset ), '%d' );
inres = inres(1);

if (seflg)
  tmp  = '<ParamLong."lBaseResolution">  {';
else
  tmp  = '<ParamLong."BaseResolution">  {';
end;
indx = strfind(asc,tmp);
outres = sscanf( asc( indx(1)+length(tmp)+(0:20) ), '%d' );
outres = outres(1);

if (seflg)
  tmp  = '<ParamDouble."aflRegridADCDuration">';
else
  tmp  = '<ParamDouble."ADCDuration">  {';
end;
indx = strfind(asc,tmp);
tmp2 = asc( indx(1)+length(tmp)+(0:100)+offset);  
indx = regexp( tmp2, '\d*' );
t_sw = sscanf( tmp2( indx(2)+(0:10) ), '%d' );

fprintf('inres: %d, outres: %d, t_ru: %d, t_ds: %d, t_ft: %d, t_sw: %d\n', ...
        inres, outres, t_ru, t_ds, t_ft, t_sw );
vrgf = mk_vrgf_dat( inres, outres, t_ru, t_ds, t_ft, t_sw );

% keyboard;

return;
