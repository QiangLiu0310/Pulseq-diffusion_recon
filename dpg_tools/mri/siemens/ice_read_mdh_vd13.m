function sMDH=ice_read_mdh_vd13(fp)

%% based on sMDH definition in n4/pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h

%
% written by W. Scott Hoge, Brigham and Women's Hosp, Boston, MA.
% (shoge -at- ieee -dot- org)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sMDH.ulDMALength=fread(fp,1,'uint32');
sMDH.lMeasUID=fread(fp,1,'uint32');
sMDH.ulScanCounter=fread(fp,1,'uint32');
sMDH.ulTimeStamp=fread(fp,1,'uint32');
sMDH.ulPMUTimeStamp=fread(fp,1,'uint32');
sMDH.ushSystemType = fread(fp,1,'uint16');
sMDH.ulPTABPosDelay = fread(fp,1,'uint16');

sMDH.ulPTABPosX=fread(fp,1,'int32');
sMDH.ulPTABPosY=fread(fp,1,'int32');
sMDH.ulPTABPosZ=fread(fp,1,'int32');

tmp = fread(fp,1,'uint32');             % reserved

sMDH.aulEvalInfoMask=fread(fp,2,'uint32');    % 8-byte, 2 elements
sMDH.ushSamplesInScan=fread(fp,1,'uint16');
sMDH.ushUsedChannels=fread(fp,1,'uint16');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sMDH.sLC.ushLine=fread(fp,1,'ushort');                %28-byte loop counter
sMDH.sLC.ushAcquisition=fread(fp,1,'ushort');
sMDH.sLC.ushSlice=fread(fp,1,'ushort');
sMDH.sLC.ushPartition=fread(fp,1,'ushort');
sMDH.sLC.ushEcho=fread(fp,1,'ushort');
sMDH.sLC.ushPhase=fread(fp,1,'ushort');
sMDH.sLC.ushRepetition=fread(fp,1,'ushort');
sMDH.sLC.ushSet=fread(fp,1,'ushort');
sMDH.sLC.ushSeg=fread(fp,1,'ushort');
sMDH.sLC.ushIda=fread(fp,1,'ushort');
sMDH.sLC.ushIdb=fread(fp,1,'ushort');
sMDH.sLC.ushIdc=fread(fp,1,'ushort');
sMDH.sLC.ushIdd=fread(fp,1,'ushort');
sMDH.sLC.ushIde=fread(fp,1,'ushort');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sMDH.sCutOff.ushPre=fread(fp,1,'uint16');
sMDH.sCutOff.ushPost=fread(fp,1,'uint16');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sMDH.ushKSpaceCentreColumn=fread(fp,1,'uint16');
sMDH.ushCoilSelect=fread(fp,1,'uint16');
% sMDH.ushDummy=fread(fp,1,'ushort');
sMDH.fReadOutOffcentre=fread(fp,1,'float');
sMDH.ulTimeSinceLastRF=fread(fp,1,'uint32');
sMDH.ushKSpaceCentreLineNo=fread(fp,1,'uint16');
sMDH.ushKSpaceCentrePartitionNo=fread(fp,1,'uint16');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slice Data:
sMDH.sSD.sSlicePosVec.flSag=fread(fp,1,'float');
sMDH.sSD.sSlicePosVec.flCor=fread(fp,1,'float');
sMDH.sSD.sSlicePosVec.flTra=fread(fp,1,'float');
sMDH.sSD.aflQuaternion=fread(fp,4,'float');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sMDH.aushIceProgramPara=fread(fp,24,'uint16');
sMDH.aushReservedPara=fread(fp,4,'uint16'); % 'ushort');

sMDH.ushApplicationCounter = fread(fp,1,'uint16');
sMDH.ushApplicationMask = fread(fp,1,'uint16');
sMDH.ulCRC = fread(fp,1,'uint32');
