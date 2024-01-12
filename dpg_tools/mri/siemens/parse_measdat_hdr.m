function hdr = parse_measdat_hdr( filename )
% updated 1/11/2021

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
      if (lngth      == 0 ) % < 1280)         %    % reached the end of a scan
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
  hdrLength = tmp;
  datStart = tmp;
end;
hdr.data_version  = data_version ;

bytes_to_skip_hdr = datStart;
asc = fread(fid,hdrLength-4,'char');
asc = char(asc');

% determine the acceleration factor:
if strcmp(data_version,'VD');
  tmp = '<ParamLong."lAccelFactPE">  {';  %%%%%  Yang,1/11/2021
  indx = strfind(asc,tmp);
  if ~isempty(indx)
    hdr.R = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );
%   else
%     hdr.R = 1;
  end;
else
  tmp = '<ParamLong."lAccelFactPE">  {';
  indx = strfind(asc,tmp);
  if ~isempty(indx)
    hdr.R = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );
  else
    hdr.R = 1;
  end;
end;

tmp = '<ParamLong."NAFLin">  {';
indx = strfind(asc,tmp);
if ~isempty(indx)
  hdr.R = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );
end;

% determine the number of acquired slices
tmp = '<ParamLong."NSlc">  {';
indx = strfind(asc,tmp);
hdr.NSlc = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );

% determine the number of phase encode lines in the ACS data
tmp = 'sPat.lRefLinesPE	 = 	';
indx = strfind(asc,tmp);
if ~isempty(indx)
    hdr.NRefLin = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );
end;

tmp = '<ParamLong."NLin">  {';
indx = strfind(asc,tmp);
if ~isempty(indx)
  hdr.NLin = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );
end;

tmp = '<ParamLong."TrueLines">  {';
indx = strfind(asc,tmp);
if isempty(indx)
  tmp = '<ParamLong."NTrueLin">  {';
  indx = strfind(asc,tmp);
end;
if ~isempty(indx)
hdr.NLin = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );
end;


tmp = '<ParamLong."NRepMeas">  {';
indx = strfind(asc,tmp);
hdr.NRep = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );

      
tmp = '<ParamLong."iMaxNoOfRxChannels">';
indx = strfind(asc,tmp);
hdr.ncoils = sscanf( asc( indx(1)+length(tmp)+10+(0:10) ), '%d' );
  hdr.FOVshift = sscanf( asc( indx(1)+length(tmp)+(3:6) ), '%d' );
tmp = '<ParamLong."relSliceNumber">  {';
indx = strfind(asc,tmp);
indx = indx(1) + +length(tmp);
tmp = sscanf( asc( indx + (0:512) ), '%d' );
hdr.slcindx = tmp( find( tmp >= 0 ) ) + 1;

tmp = 'sSliceAcceleration.lMultiBandFactor';
indx = strfind(asc,tmp);
if ~isempty(indx)
  hdr.SMS = sscanf( asc( indx(end)+length(tmp)+(3:6) ), '%d' );
end;
tmp = 'sSliceAcceleration.lFOVShiftFactor';
indx = strfind(asc,tmp);
if ~isempty(indx)
  hdr.FOVshift = sscanf( asc( indx(1)+length(tmp)+(3:6) ), '%d' );
end;

%%% for the VB17 data:
tmp = 'sWiPMemBlock.adFree[1]';
indx = strfind(asc,tmp);
if ~isempty(indx)
  hdr.SMS = sscanf( asc( indx(end)+length(tmp)+17+(3:6) ), '%d' );
end;
tmp = 'sWiPMemBlock.adFree[2]';
indx = strfind(asc,tmp);
if ~isempty(indx)
  hdr.FOVshift = sscanf( asc( indx(1)+length(tmp)+17+(3:6) ), '%d' );
end;


tmp = 'WiPMemBlock.alFree[4]                   =';
indx = strfind(asc,tmp);
if ~isempty(indx)
indx = indx(1) + +length(tmp);
tmp = sscanf( asc( indx + (0:512) ), '%d' );
hdr.seg = tmp;
end;




% keyboard;
% vrgf = regrid_epi_rampsamp( asc );

if ( ~isfield( hdr, 'NRefLin') && isfield( hdr, 'NLin') )
  hdr.NRefLin = hdr.NLin;
end;

fclose(fid);

return;

function vrgf = regrid_epi_rampsamp( asc )

vrgf = 1;

tmp  = '<ParamLong."FlattopTime">  {';
indx = strfind(asc,tmp);
if isempty(indx), return; end;
t_ft = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );

tmp  = '<ParamLong."RampupTime">  {';
indx = strfind(asc,tmp);
t_ru = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );

tmp  = '<ParamLong."RampdownTime">  {';
indx = strfind(asc,tmp);
t_rd = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );

tmp  = '<ParamLong."DelaySamplesTime">  {';
indx = strfind(asc,tmp);
t_ds = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );

tmp  = '<ParamLong."DestSamples">  {';
indx = strfind(asc,tmp);
inres = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );

tmp  = '<ParamLong."BaseResolution">  {';
indx = strfind(asc,tmp);
outres = sscanf( asc( indx(1)+length(tmp)+(0:10) ), '%d' );

tmp  = '<ParamDouble."ADCDuration">  {';
indx = strfind(asc,tmp);
tmp2 = asc( indx(1)+length(tmp)+(0:100));  
indx = regexp( tmp2, '\d*' );
t_sw = sscanf( tmp2( indx(2)+(0:10) ), '%d' );

vrgf = mk_vrgf_dat( inres, outres, t_ru, t_ds, t_ft, t_sw );

% keyboard;

return;
