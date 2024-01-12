function varargout = read_meas_prot(wildfile, varargin)
%READ_MEAS_PROT  read protocol from VB- and VD-style "meas.dat"
%
% YAPS = read_meas_prot(filename)


% Copyright © 2006-2012 Jonathan R. Polimeni and
%   The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense


% jonathan polimeni <jonp@padkeemao.nmr.mgh.harvard.edu>, 2007/jan/03
% $Id: read_meas_prot.m,v 1.17 2014/02/17 21:53:56 jonp Exp $
%**************************************************************************%

  VERSION = '$Revision: 1.17 $';
  if ( nargin == 0 ), help(mfilename); return; end;

  header = '';
  if ( nargin >= 2 ),
    if ( ~isstruct(varargin{1}) ),
      header = varargin{1};
    end;
  end;

  do_oldway = 0;
  ReturnDependencyMeas = 0;
  if ( nargin >= 3 ),
    ReturnDependencyMeas = varargin{2};
  end;
  DO__RETURN_DEPENDENCY_MEAS = ReturnDependencyMeas;

  IS__MULTIRAID = 0;
  
  
  measfile = mrir_sysutil__wildfile(wildfile);
  [pathstr, filestr, extstr] = fileparts(measfile);


  %------------------------------------------------------------------------%
  % error checking

  if ( ~exist(measfile, 'file') ),
    error('file [%s] does not exist', measfile);
  end;

  VERBOSE = 0;
  IS__VD_PLATFORM = 0;


  %------------------------------------------------------------------------%

  if ( isempty(header) ),
    [fp, errstr] = fopen(measfile, 'r', 'l');
    if ( fp == -1 ),
      error(errstr);
    end;

    % determine size (in bytes) of ascii header files stored in the 'meas.dat'
    % format (i.e., "Config_.evp", "Dicom_.evp", etc) to skip over them all.
    % [note: this works for VB11A also---the first integer is 32, which is
    % the number of bytes to skip at the beginning!]
    data_start = fread(fp, 1, 'uint32');

    if ( data_start == 0 ),    % Is this a VD/VE data file?
      IS__VD_PLATFORM = 1;
    end;

    if ( do_oldway ),

      if ( IS__VD_PLATFORM ),                                % Read the real data offset
        fseek(fp, 16, 'bof');
        header_start = fread(fp, 1, 'uint32');
        fseek(fp, header_start, 'bof');
        data_start = header_start + fread(fp, 1, 'uint32');
      end;

      header_files = fread(fp, 1, 'uint32');

      % read header into one string for parsing
      header = fscanf(fp, '%c', data_start-24);

      % HACK: sometimes two headers are stored in one file (e.g.,
      % AdjCoilSens); find the end the YAPS marker and if there are two
      % truncate header at the first one
      header_end = regexp(header, '### ASCCONV BEGIN');  % avoids MeasYaps
      if ( length(header_end) > 1 ),
        header = header(1:header_end(end));
      end;

    else,  % newway

      if ( IS__VD_PLATFORM ),

        % count number of concatenated meas.dat files stored in this file
        num_dependent_files = fread(fp, 1, 'uint32');

        if ( num_dependent_files > 1 ),
          IS__MULTIRAID = 1;
          fprintf(1, '\n ***MULTI-RAID file detected***\n  %d meas contained within file \n\n', num_dependent_files);

          if ( DO__RETURN_DEPENDENCY_MEAS == 0 ),
            fprintf(1, ' input "%s" contains [%d] dependent meas.dat ADJUSTMENT file(s)\n', strcat(filestr, extstr), num_dependent_files-1);
            disp(' reading in *main data file* in multi-RAID file...');
            if ( num_dependent_files == 2 ),
              fprintf(1, ' re-run and set input argument "ReturnDependencyMeas" to 1 to load another meas within this file\n\n');
            elseif ( num_dependent_files == 3 ),
              fprintf(1, ' re-run and set input argument "ReturnDependencyMeas" to 1 or 2 to load another meas within this file\n\n');
            else,
              fprintf(1, ' re-run and set input argument "ReturnDependencyMeas" to 1--%d to load another meas within this file\n\n', num_dependent_files-1);
            end;
          else,
            fprintf(1, ' reading prot %d of %d within the multi-RAID file\n\n', ReturnDependencyMeas, num_dependent_files);
          end;
        end;

        if ( ~IS__MULTIRAID ),
          meas_file_within_multifile = 0;
        end;

        % by default read the last file, which is the main file stored after the adjustment files
        if ( DO__RETURN_DEPENDENCY_MEAS == 0 ),
          meas_file_within_multifile = (num_dependent_files - 1);
        else,
          meas_file_within_multifile = ReturnDependencyMeas - 1;
        end;

        % magic number alert: 152
        fseek(fp, 16+(meas_file_within_multifile*152), 'bof');
        header_start = fread(fp, 1, 'uint64');
        fseek(fp, header_start, 'bof');
        data_start = header_start + fread(fp, 1, 'uint32');

        % read header into one string for parsing
        header = char(fread(fp, data_start-header_start-4, 'uchar').');

      else,
        header_files = fread(fp, 1, 'uint32');

        % read header into one string for parsing
        header = fscanf(fp, '%c', data_start-24);

        % HACK: sometimes two headers are stored in one file (e.g.,
        % AdjCoilSens); find the end the YAPS marker and if there are two
        % truncate header at the first one
        header_end = regexp(header, '### ASCCONV BEGIN');  % avoids MeasYaps
        if ( length(header_end) > 1 ),
          header = header(1:header_end(end));
        end;

      end;  %% IS__VD_PLATFORM


    end;

    fclose(fp);
  end;   %% HEADER


  %------------------------------------------------------------------------%
  % ICE parameters from Config_.evp

  param_list = {
      ... %%% stolen from MDH (useful for asymmetric echo acquisitions)
      'MDH_ushSamplesInScan', ...
      ... %%% k-space sample counts  (no "RawSeg" for some reason)
      'RawCol','RawLin','RawCha','RawSet','RawEco','RawPhs','RawRep', ...
      'RawPar','RawSlc','RawIda','RawIdb','RawIdc','RawIdd','RawIde','RawAve', ...
      ... %%%
      ... %%% iPAT / parallel imaging acquisition parameters
      'NAFLin','NAFPar','NFirstLin','NFirstPar', ...
      'NRefLin','NRefPar','NFirstRefLin','NFirstRefPar', ...
      ... %%%
      ... %%% effective k-space dimensions
      'NColMeas','NLinMeas','NChaMeas','NSetMeas','NEcoMeas','NPhsMeas','NRepMeas','NSegMeas', ...
      'NParMeas','NSlcMeas','NIdaMeas','NIdbMeas','NIdcMeas','NIddMeas','NIdeMeas','NAveMeas', ...
      ... %%%
      ... %%% final image dimensions
      'NImageCols','NImageLins','NImagePar','NImagePars'};

  Config_evp = cell2struct(cell(length(param_list),1), param_list, 1);
  Config_evp_values = [];

  for ind = 1:length(param_list),
    param = param_list{ind};


    %%%%%%% MATCH TYPE #0: integer dimension size, e.g., <ParamLong."RawCol">  { 256  }
    match = regexp(header, ['<Param\w+\."(?<param>' param ')">\s*{\s*(?<value>[-]*\d*)\s*}'], 'names');

    % check if no match is found (can happen naturally if, e.g., no accel no NRefLin found)
    if ( isempty(match) ),
      if ( VERBOSE ),
        warning('SIEMENS:IO:versioning', 'missing protocol parameter "%s"---check sequence', param);
      end;
      continue;
    end;

    % [[jrp, 2011/sep/02]] consider second or first match; was always first!
    % (workaround for MEMPRAGE & partition/inner GRAPPA w/moco dev)
    if ( length(match) > 1 ),
      match = match(end);
    else,
      match = match(1);
    end;

    % empty means value for this parameter dimension = 0
    if ( isempty(match.value) ),
      match.value = '0';
    end;

    % save out struct and numerical array
    param_values(ind) = str2num(match.value);

    %%%-------------------------------%
    %%% SPECIAL CASES: "NFirstLin", "NFirstPar", "NFirstRefLin", "NFirstRefPar" (0-based indices needing conversion to 1-based)
    if ( regexp(param, '^NFirst') ),
      param_values(ind) = param_values(ind) + 1;
    end;


    Config_evp = setfield(Config_evp, param_list{ind}, param_values(ind));

  end;  %% FOR

  %  %%%-------------------------------%
  %  %%% SPECIAL CASES: "NFirstRefLin", "NFirstRefPar"
  %  if ( Config_evp.NFirstRefLin == 0 ),
  %    Config_evp.NRefLin = 0;
  %  end;
  %  if ( Config_evp.NFirstRefPar == 0 ),
  %    Config_evp.NRefPar = 0;
  %  end;


  %========================================================================%
  %------------------------------------------------------------------------%
  %========================================================================%
  % IDEA parameters from YAPS, a.k.a., ASCCONV (contained in meas.asc)

  %  TODO: determine phase encoding or readout direction, and polarities,
  %  to convert reconstructed matrix into proper viewing position

  %------------------------------------------------------------------------%
  % anastasia's list

  %  p.SequenceFileName      = 'tSequenceFileName';
  %  p.ProtocolName          = 'tProtocolName';
  %  p.FlipAngleDegree       = 'adFlipAngleDegree';
  %  p.TxFrequency           = 'sTXSPEC.lFrequency';
  %  p.BaseResolution        = 'sKSpace.lBaseResolution';
  %  p.PhaseEncodingLines    = 'sKSpace.lPhaseEncodingLines';
  %  p.FourierRows           = 'iNoOfFourierLines';
  %  p.FourierColumns        = 'iNoOfFourierColumns';
  %  p.NumChannels           = 'iMaxNoOfRxChannels';
  %  p.PhaseResolution       = 'sKSpace.dPhaseResolution';
  %  p.SliceResolution       = 'sKSpace.dSliceResolution';
  %  p.Partitions            = 'sKSpace.lPartitions';
  %  p.TR                    = 'alTR\[[0-9]*\]';
  %  p.TI                    = 'alTI\[[0-9]*\]';
  %  p.TE                    = 'alTE\[[0-9]*\]';
  %  p.SliceArraySize        = 'sSliceArray.lSize';
  %  p.SlicePositionSag      = 'sSliceArray.asSlice\[[0-9]*\].sPosition.dSag';
  %  p.SlicePositionCor      = 'sSliceArray.asSlice\[[0-9]*\].sPosition.dCor';
  %  p.SlicePositionTra      = 'sSliceArray.asSlice\[[0-9]*\].sPosition.dTra';
  %  p.SliceNormalSag        = 'sSliceArray.asSlice\[[0-9]*\].sNormal.dSag';
  %  p.SliceNormalCor        = 'sSliceArray.asSlice\[[0-9]*\].sNormal.dCor';
  %  p.SliceNormalTra        = 'sSliceArray.asSlice\[[0-9]*\].sNormal.dTra';
  %  p.Thickness             = 'sSliceArray.asSlice\[[0-9]*\].dThickness';
  %  p.PhaseFOV              = 'sSliceArray.asSlice\[[0-9]*\].dPhaseFOV';
  %  p.ReadoutFOV            = 'sSliceArray.asSlice\[[0-9]*\].dReadoutFOV';
  %  p.InPlaneRot            = 'sSliceArray.asSlice\[[0-9]*\].dInPlaneRot';
  %  p.DiffDirections        = 'sDiffusion.lDiffDirections';
  %  p.bValue                = 'sDiffusion.alBValue\[1\]';
  %  p.b0Repetitions         = 'sWiPMemBlock.alFree\[8\]';
  %  p.Repetitions           = 'lRepetitions';


  % dwelltime = The base resolution and the effective dwelltime are specified without oversampling
  % (true number of samples = base_resolution * oversampling,
  % real dwelltime = effective dwelltime/oversampling).

  % duration of the ADC (columns * dwelltime)


  % YAPS.lRampTime, YAPS.lFlatTime, YAPS.lADCDuration, YAPS.lRampMode

  [YAPS, param_list] = read_meas_prot__struct;
  YAPS.file = measfile;


  param_values = zeros(length(param_list),1);


  % scan through header for each of the protocol parameter values
  for ind = 1:length(param_list),
    param = param_list{ind};


    %%%-------------------------------%
    %%% SPECIAL CASE: coil element names -- returns an array of coil name strings

    if ( strcmp(param, 'sCoilElementID_tElement') ),

      match = regexp(header, ['asCoilSelectMeas\S*(?<param>lRxChannelConnected)\s*=\s*(?<value>\d*\.?x?\d*)'], 'names');
      % note this is tricksy: VD-platform uses (for some reason)
      % "sCoilSelectMeas" at beginning of string whereas in prior versions
      % "asCoilSelectMeas" is used! so the above will not match VD-platform
      % files.

      if ( isempty(match) ),
        % VD-platform uses ADCs rather than RF channels since two coils
        % share one receiver line
        match = regexp(header, ['sCoilSelectMeas\S*(?<param>lADCChannelConnected)\s*=\s*(?<value>\d*\.?x?\d*)'], 'names');
      end;

      channel_val = str2double({match.value});
      [channel_unique, channel_index] = unique(channel_val);


      % need the "Rx" in there since Tx coils are also listed (VD-platform!)
      match = regexp(header, ['sCoilSelectMeas\.aRx\S*(?<param>tElement)\s*=\s*"(?<string>\S*)"'], 'names');

      if ( isempty(match) ),
        % (pre-VD-platform)
        match = regexp(header, ['sCoilSelectMeas\S*(?<param>tElement)\s*=\s*"(?<string>\S*)"'], 'names');
      end;

      % sorted by the receiver channel ID number (not alphabetically by
      % element string), so should match native ordering of image data.
      % NOTE that order listed in "sCoilSelectMeas" in header does not
      % necessarily match order found in MDHs. need to sort by channel ID
      % number
      YAPS = setfield(YAPS, 'sCoilElementID_tElement', {match(channel_index).string});
      continue;

    end;


    %%%-------------------------------%
    %%% SPECIAL CASE: tCoilID

    % sCoilSelectMeas.aRxCoilSelectData[0].asList[0].sCoilElementID.tCoilID    =      "Brain_64"

    if ( strcmp(param, 'tCoilID') ),

      match = regexp(header, ['asCoilSelectMeas\S*(?<param>tCoilID)\s*=\s*"(?<string>\S*)"'], 'names');

      if ( ~isempty(match) ),
        % consider only first match
        match = match(1);

        % set string value
        YAPS = setfield(YAPS, regexprep(param_list{ind}, '\.', '_'), match.string);
        continue;
      end;

    end;


    %%%-------------------------------%
    %%% SPECIAL CASE: slice array params, "dPhaseFOV", "dReadoutFOV"

    if ( strcmp(param, 'dPhaseFOV') || strcmp(param, 'dReadoutFOV') ),
      param = strcat('asSlice\S*', param);
    end;


    %%%-------------------------------%

    % exploit MATLAB regexp machinery to pull out parameter/value pairs

    %%%%%%% MATCH TYPE #1: string parameter, e.g.,  tSequenceFileName                        = "%SiemensSeq%\gre"
    match = regexp(header, ['\W(?<param>' param ')\s*=\s*"(?<string>\S*)"'], 'names');
    if ( ~isempty(match) ),
      % consider only first match
      match = match(1);

      % set string value
      YAPS = setfield(YAPS, regexprep(param_list{ind}, '\.', '_'), match.string);
      continue;
    end;



    %%%%%%% MATCH TYPE #2a: numerical parameter but single element, e.g., iMaxNoOfRxChannels                       = 9
    match = regexp(header, ['\W(?<param>' param ')\s*=\s*(?<value>(-*)\d*\.?x?\w*)'], 'names', 'once');

    if ( isempty(match) ),
      %%%%%%% MATCH TYPE #2b: single-element numerical parameter, but with <Precision> tag
      match = regexp(header, ['<Param\w*\."(?<param>' param ')">\s*\{\s*<Precision> \d+ \s(?<value>(-*)\d*\.?x?\w*)'], 'names', 'once');
    end;


    if ( isempty(match) ),


      %%%%%%% MATCH TYPE #3: numerical parameter array, possibly negative, e.g., sRXSPEC.alDwellTime[0]                   = 6500
      match = regexp(header, ['\W(?<param>' param ')\[(?<index>\d*)\]\s*=\s*(?<value>(-*)\d*\.?\d*)'], 'names');

      if ( ~isempty(match) ),

        %%%-------------------------------%
        %%% SPECIAL CASE: "alRegridMode"

        % the "alRegridMode" parameter consists of an array of values indicating
        % which basic waveform components comprise the readout waveform. the
        % presence of a fifth component (i.e., "alRegridMode[4]") indicates a
        % sinusoidal waveform. the parameter value itself is irrelevant here.
        if ( strcmp(param, 'alRegridMode') ),
          sinusoid_gradient_waveform_flag = 0;
          for im = 1:length(match),
            if ( match(im).index == 4 ),
              % sinusoidal waveform found!
              % save out struct and numerical array
              sinusoid_gradient_waveform_flag = 1;
            end;
          end;
          if ( sinusoid_gradient_waveform_flag == 1 ),
            param_values(ind) = 4;  % 'REGRID_SINUSOIDAL'
          else,
            param_values(ind) = 2;  % 'REGRID_TRAPEZOIDAL'
          end;
          YAPS = setfield(YAPS, param_list{ind}, param_values(ind));
          continue;
        end;

        %%%-------------------------------%


        indices = str2num([match.index].');
        values = [];

        for ii = 1:length(match),
          % empty means value for this parameter = 0
          if ( isempty(match(ii).value) ),
            match(ii).value = '0';
          end;
          values(str2num(match(ii).index)+1) = str2num(match(ii).value);
        end;

        YAPS = setfield(YAPS, regexprep(param_list{ind}, '\.', '_'), values);
        continue;

      end;

    end;

    if ( isempty(match) ),
      %%%%%%% MATCH TYPE #5: find a string/value pair in the mark-up---this may be all there is in the future!
      match = regexp(header, ['<ParamString\."(?<param>' param ')">\s*{\s*"(?<string>[^\t\r\n\f\v]*)"'], 'names');

      if ( ~isempty(match) ),
        % consider only first match
        match = match(1);

        % set string value
        YAPS = setfield(YAPS, regexprep(param, '\.', '_'), regexprep(match.string, ' ', '_'));
        continue;
      end;
    end;

    if ( isempty(match) ),
      %%%%%%% MATCH TYPE #6: find a long/value pair in the mark-up---this may be all there is in the future!
      match = regexp(header, ['<ParamLong\."(?<param>' param ')">\s*{\s*(?<value>\d*)\s*}'], 'names');

      if ( ~isempty(match) ),
        % consider only first match
        match = match(1);

        % set string value
        YAPS = setfield(YAPS, regexprep(param, '\.', '_'), str2num(match.value));
        continue;
      end;
    end;

    % check if still no match is found
    if ( isempty(match) ),
      if ( VERBOSE ),
        warning('SIEMENS:IO:versioning', 'missing protocol parameter "%s"---check sequence', param);
      end;

      if ( regexp(param, '^fl.*') ),
        YAPS = setfield(YAPS, regexprep(param_list{ind}, '\.', '_'), 0.0);
      end;

      continue;
    end;


    %========================================================================%
    %========================================================================%
    % consider only first match
    match = match(1);

    % empty means value for this parameter = 0
    if ( isempty(match.value) ),
      match.value = '0';
    end;

    % convert hex entries to decimal
    if ( regexp(match.value, '^0x\d*') ),
      match.value = num2str(hex2dec( regexprep(match.value, '^0x', '') ));
    end;


    %%%-------------------------------%
    %%% ENUMERATED TYPE: "ucPhasePartialFourier" and "ucSlicePartialFourier"
    if ( strcmp(param, 'ucPhasePartialFourier') || strcmp(param, 'ucSlicePartialFourier') ),

      % "enum PartialFourierFactor" snarfed from <n4/pkg/MrServers/MrProtSrv/MrProt/SeqDefines.h>
      switch str2num(match.value),
       case hex2dec('01'),  % PF_HALF = 0x01
        match.value = '4/8';
       case hex2dec('02'),  % PF_5_8  = 0x02
        match.value = '5/8';
       case hex2dec('04'),  % PF_6_8  = 0x04
        match.value = '6/8';
       case hex2dec('08'),  % PF_7_8  = 0x08
        match.value = '7/8';
       case hex2dec('10'),  % PF_OFF  = 0x10
        match.value = '8/8';
       otherwise,
        warning(sprintf('unrecognized mode for "%s":  [ %s ]', param, match.value));
      end;

    end;


    %%%-------------------------------%
    %%% ENUMERATED TYPE: "ucMultiSliceMode"

    if ( strcmp(param, 'ucMultiSliceMode') ),
      % "enum MultiSliceMode" snarfed from <n4/pkg/MrServers/MrProtSrv/MrProt/SeqDefines.h>
      switch str2num(match.value),
       case hex2dec('01'),  % MSM_SEQUENTIAL  = 0x01
        match.value = 'MSM_SEQUENTIAL';
       case hex2dec('02'),  % MSM_INTERLEAVED = 0x02
        match.value = 'MSM_INTERLEAVED';
       case hex2dec('04'),  % MSM_SINGLESHOT  = 0x04
        match.value = 'MSM_SINGLESHOT';
       otherwise,
        warning(sprintf('unrecognized mode for "%s":  [ %s ]', param, match.value));
      end;

      YAPS = setfield(YAPS, regexprep(param_list{ind}, '\.', '_'), match.value);
      continue;
    end;


    %%%-------------------------------%
    %%% ENUMERATED TYPE: "ucSliceSeriesMode"

    if ( strcmp(param, 'sSliceArray.ucMode') ),
      % "enum ???" snarfed from <n4/pkg/MrServers/MrProtSrv/MrProt/SeqDefines.h>
      switch str2num(match.value),
       case hex2dec('01'),  % ASCENDING   = 0x01
        match.value = 'ASCENDING';
       case hex2dec('02'),  % DESCENDING  = 0x02
        match.value = 'DESCENDING';
       case hex2dec('04'),  % INTERLEAVED = 0x04
        match.value = 'INTERLEAVED';
       case hex2dec('08'),  % AUTOMATIC   = 0x08
        match.value = 'AUTOMATIC';
       case hex2dec('10'),  % APEXTOBASE  = 0x10
        match.value = 'APEXTOBASE';
       case hex2dec('20'),  % BASETOAPEX  = 0x20
        match.value = 'BASETOAPEX';
       otherwise,
        warning(sprintf('unrecognized mode for "%s":  [ %s ]', param, match.value));
      end;

      YAPS = setfield(YAPS, regexprep(param_list{ind}, '\.', '_'), match.value);
      continue;
    end;


    %%%-------------------------------%
    %%% ENUMERATED TYPE: "ucPATMode"
    if ( strcmp(param, 'ucPATMode') ),

      % "enum PATSelMode" snarfed from <n4/pkg/MrServers/MrProtSrv/MrProt/SeqDefines.h>
      switch str2num(match.value),
       case hex2dec('01'),
        match.value = 'PAT_MODE_NONE';
       case hex2dec('02'),
        match.value = 'PAT_MODE_GRAPPA';
       case hex2dec('04'),
        match.value = 'PAT_MODE_SENSE';
       case hex2dec('08'),
        match.value = 'PAT_MODE_2D';
       case hex2dec('10'),
        match.value = 'PAT_MODE_CAIPI';
       otherwise,
        warning(sprintf('unrecognized mode for "%s":  [ %s ]', param, match.value));
      end;
      
%      <ParamChoice."ucPATMode"> 
%        {
%          <Label> "iPAT Mode" 
%          <Tooltip> "Modes for iPAT" 
%          <Visible> "true" 
%          <Default> "GRAPPA" 
%          <Limit> { "GRAPPA" "mSense" "CAIPIRINHA" }
%        }

%        <ParamLong."ulCaipirinhaMode"> 
%        {
%          <Default> 1 
%          <Limit> { 1 2 4 8 16 }
%        }
%        
%        <ParamLong."lReorderingShift3D"> 
%        {
%        }



      YAPS = setfield(YAPS, regexprep(param_list{ind}, '\.', '_'), match.value);
      continue;
    end;


    %%%-------------------------------%
    %%% ENUMERATED TYPE: "ucRefScanMode"
    if ( strcmp(param, 'ucRefScanMode') ),

      % "enum PATRefScanMode" snarfed from <n4/pkg/MrServers/MrProtSrv/MrProt/SeqDefines.h>
      switch str2num(match.value),
       case hex2dec('01'),
        match.value = 'PAT_REF_SCAN_UNDEFINED';      % e.g. if no PAT is selected
       case hex2dec('02'),
        match.value = 'PAT_REF_SCAN_INPLACE';        % sequence supplies inplace reference lines
       case hex2dec('04'),
        match.value = 'PAT_REF_SCAN_EXTRA';          % sequence supplies extra reference lines
       case hex2dec('08'),
        match.value = 'PAT_REF_SCAN_PRESCAN';        % sequence does not supply reference lines, the data must have been acquired with a previous measurement
       case hex2dec('10'),
        match.value = 'PAT_REF_SCAN_INTRINSIC_AVE';  % The sequence contains intrinsic ref.lines due to sharing e.g. in the averages dimension
       case hex2dec('20'),
        match.value = 'PAT_REF_SCAN_INTRINSIC_REP';  % The sequence contains intrinsic ref.lines due to sharing e.g. in the repetition or phases dimension (i.e., TSENSE)
       case hex2dec('40'),
        match.value = 'PAT_REF_SCAN_INTRINSIC_PHS';  % The sequence contains intrinsic ref.lines due to sharing e.g. in the repetition or phases dimension (i.e., TSENSE)
       case hex2dec('80'),
        match.value = 'PAT_REF_SCAN_INPLACE_LET';    % A single (L)ong (E)cho (T)rain acquires reference lines and imaging lines
       otherwise,
        warning(sprintf('unrecognized mode for "%s":  [ %s ]', param, match.value));
      end;

      YAPS = setfield(YAPS, regexprep(param_list{ind}, '\.', '_'), match.value);
      continue;
    end;


    %------------------------------------------------------------------------%

    % save out struct and numerical array
    param_values(ind) = str2num(match.value);
    YAPS = setfield(YAPS, regexprep(param_list{ind}, '\.', '_'), param_values(ind));

  end;  %% FOR

  %------------------------------------------------------------------------%
  %------------------------------------------------------------------------%
  %%% post-processing


  if ( ~strncmp(YAPS.ulVersion, '0x', 2) ),
    YAPS.ulVersion = lower(['0x' dec2hex(YAPS.ulVersion)]);
  end;

  switch YAPS.ulVersion,
   case '0xbee332',
    YAPS.ulVersion = [YAPS.ulVersion, ':  "N4_VA25A_LATEST_20040724"'];
   case '0x1421cf5',
    YAPS.ulVersion = [YAPS.ulVersion, ':  "N4_VB11D_LATEST_20040724"'];
   case '0x1452a3b',
    YAPS.ulVersion = [YAPS.ulVersion, ':  "N4_VB13A_LATEST_20060607"'];
   case '0x1483779',
    YAPS.ulVersion = [YAPS.ulVersion, ':  "N4_VB15A_LATEST_20070519"'];
   case '0x14b44b6',
    YAPS.ulVersion = [YAPS.ulVersion, ':  "N4_VB17A_LATEST_20090307"'];
   case '0x273bf24',
    YAPS.ulVersion = [YAPS.ulVersion, ':  "N4_VD11D_LATEST_20110129"'];
   case '0x276a554',
    YAPS.ulVersion = [YAPS.ulVersion, ':  "N4_VD13C_LATEST_20121124"'];
   case '0x30c0783',
    YAPS.ulVersion = [YAPS.ulVersion, ':  "N4_VE11B_LATEST_20150530"'];
   case '0x30c2e91',
    YAPS.ulVersion = [YAPS.ulVersion, ':  "N4_VE11C_LATEST_20160120"'];
   otherwise,
    disp(sprintf('==> [%s]: version string "%s" not on record!', measfile, ...
                 YAPS.ulVersion));
  end;


  %%%-------------------------------%
  %%% SPECIAL CASE: "alTE"

  YAPS.alTE = YAPS.alTE(1:YAPS.lContrasts);


  %%%-------------------------------%
  %%% SPECIAL CASE: regrid parameters (only a pain in VB15)

  if ( isempty(YAPS.alRegridMode) ),

    param_list = {'alRegridRampupTime', 'alRegridRampdownTime', ...
                  'alRegridFlattopTime', 'alRegridDelaySamplesTime', ...
                  'alRegridMode', 'aflRegridADCDuration'};

    for ind = 1:length(param_list),
      param = param_list{ind};

      match = regexp(header, ['<Param\w+\."(?<param>' param ')">\s*{\s*(?<value>[\d ]*)\s*}'], 'names');

      % SPECIAL CASE (within a special case): pesky 'Precision' string for "aflRegridADCDuration"
      if ( isempty(match) ),
        match = regexp(header, ['<Param\w+\."(?<param>' param ')">\s*{\s*<Precision>[ \d]*\s*(?<value>[\d \.]*)\s*}'], 'names');
      end;

      match.value = str2num(match.value);
      YAPS = setfield(YAPS, param_list{ind}, match.value(1));

    end;

  end;


  %%%-------------------------------%
  %%% ENUMERATED TYPE: "alRegridMode"

  switch YAPS.alRegridMode,
    % "enum RegriddingMode" snarfed from <n4/pkg/MrServers/MrProtSrv/MrProt/SeqDefines.h>

   case hex2dec('01'),
    YAPS.ucRegridMode = 'REGRID_NONE';
   case hex2dec('02'),
    YAPS.ucRegridMode = 'REGRID_TRAPEZOIDAL';
   case hex2dec('04'),
    YAPS.ucRegridMode = 'REGRID_SINUSOIDAL';
   otherwise,
    warning(sprintf('unrecognized mode for "%s":  [ %d ]', 'alRegridMode', YAPS.alRegridMode));
  end;


  %%%-------------------------------%
  %%% SPECIAL CASE: flBandwidthPerPixelPhaseEncode and iEffectiveEpiEchoSpacing (only present in VB15)

  param_list = {'flBandwidthPerPixelPhaseEncode', 'iEffectiveEpiEchoSpacing'};

  for ind = 1:length(param_list),
    param = param_list{ind};

    match = regexp(header, ['<Param\w+\."(?<param>' param ')">\s*{\s*(?<value>[\d.?\d]*)\s*}'], 'names');

    if ( isempty(match) ),
      continue;
    end;

    match.value = str2num(match.value);

    if ( isempty(match.value) ),
      match.value = 0;
    end;

    YAPS = setfield(YAPS, param_list{ind}, match.value(1));

  end;


  %%%-------------------------------%
  %%% SPECIAL CASE: sSliceArray

  param_list = {'sPosition.dSag', 'sPosition.dCor', 'sPosition.dTra', ...
                  'sNormal.dSag',   'sNormal.dCor',   'sNormal.dTra', ...
                'dThickness', 'dPhaseFOV', 'dReadoutFOV', 'dInPlaneRot'};

  % empty sub-struct for slice array
  asSlice = cell2struct(cell(length(param_list),1), regexprep(param_list, '\.', '_'), 1);

  for slc = 1:YAPS.sSliceArray_lSize,

    % assign new empty sub-struct for each slice
    YAPS.sSliceArray(slc) = asSlice;

    for ind = 1:length(param_list),

    param = param_list{ind};
    param = regexprep(param, '\.', '\\.');

    match = regexp(header, ['sSliceArray\.asSlice\[' num2str(slc-1) '\]\.(?<param>' param ')\s*=\s*(?<value>(-*)\d*\.?x?\d*(e-?\+?\d+)?)'], 'names');

    % check if no match is found (e.g., only one of three possible normals will be represented, so empty match is common)
    if ( isempty(match) ),
      match(1).value = '0';
    end;

    % consider only first match
    match = match(1);


    YAPS.sSliceArray(slc) = setfield(YAPS.sSliceArray(slc), regexprep(param_list{ind}, '\.', '_'), str2num(match.value));
    end;
  end;


  %%%-------------------------------%
  %%% SPECIAL CASE: slice ordering

%   % AnatomicalSliceNo
%   % ProtocolSliceNumber
%             <ParamArray."ProtocolSliceNumber">
%             {
%               <MinSize> 0
%               <MaxSize> 2147483647
%               <Default> <ParamLong."">
%               {
%               }
%               { }    <--- first entry is "0", but by convention this is left as blank
%               { 2  }
%               { 4  }
%               { 6  }
%               { 8  }
%               { 10  }
%               { 12  }
%               { 14  }
%               { 16  }
%               { 18  }
%               { 20  }
%               { 22  }
%               { 24  }
%               { 26  }
%               { 28  }
%               { 1  }
%               { 3  }
%               { 5  }
%               { 7  }
%               { 9  }
%               { 11  }
%               { 13  }
%               { 15  }
%               { 17  }
%               { 19  }
%               { 21  }
%               { 23  }
%               { 25  }
%               { 27  }
%               { -1  }
%               { -1  }
%               { -1  }
%               { -1  }
%               { -1  }
%
%
%   % relSliceNumber
%   % chronSliceIndices
%
%         <ParamLong."chronSliceIndices">
%         {
%           <Default> -1
%           0 2 4 6 8 10 12 14 16 18
%           20 22 24 26 28 1 3 5 7 9
%           11 13 15 17 19 21 23 25 27 -1
%           -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
%           -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
%           -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
%           -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
%           -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
%           -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
%           -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
%           -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
%           -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
%           -1 -1 -1 -1 -1 -1 -1 -1
%         }



  %------------------------------------------------------------------------%
  % derived fields

  % (1) "flBandwidth": siemens bandwidth, as specified in protocol
  % (which is readout bandwidth, not pixel bandwidth)

  for contrast = 1:YAPS.lContrasts,

    % sample dwell time is given in nanoseconds in header
    if ( ~isempty(Config_evp.NColMeas) ),
      bandwidth = 1 / (YAPS.alDwellTime(contrast) * Config_evp.NColMeas / 1e9);
      YAPS.flBandwidthPerPixel__DERIVED(contrast) = bandwidth;
    end;
    YAPS.flFullReadoutBandwidth__DERIVED(contrast) = 1 / (YAPS.alDwellTime(contrast) / 1e9);

  end;


  %------------------------------------------------------------------------%

  if ( nargout > 0 ),
    varargout{1} = YAPS;
    varargout{2} = Config_evp;
    varargout{3} = header;
  end;


  return;


  %************************************************************************%
  %%% $Source: /space/padkeemao/1/users/jonp/cvsjrp/PROJECTS/IMAGE_RECON/mrir_toolbox/read_meas_prot.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
