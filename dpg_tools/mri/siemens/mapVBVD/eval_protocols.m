function tab = eval_protocols(pathstr)
    
    search_str = '*.dat';
    if exist('pathstr','var')
        search_str = fullfile(pathstr, search_str);
    end
    filelist = dir(search_str);
    
    for k=1:numel(filelist)
        prot = eval_twix_hdr(fullfile(pathstr,filelist(k).name));
        
        file{k} = filelist(k).name;
        seq{k} = prot.Meas.SequenceDescription;
        mid{k} = prot.Meas.MeasUID;

        FA{k} = prot.MeasYaps.adFlipAngleDegree{1};
        voltage{k} = prot.MeasYaps.sTXSPEC.asNucleusInfo{1}.flReferenceAmplitude;
        
        TR{k} = prot.MeasYaps.alTR{1}/1000;
        TE{k} = prot.MeasYaps.alTE{1}/1000;

        TI{k} = prot.Meas.alTI(1)/1e6;
        if numel(prot.MeasYaps.alTR)>1
            TR2{k} = prot.MeasYaps.alTR{2}/1e6;
        else
            TR2{k} = [];
        end

        TRvol{k} = TR{k} * round(1e6 *prot.MeasYaps.lScanTimeSec/prot.MeasYaps.alTR{1}/prot.Meas.NSetMeas)/1000;
        rep{k} = prot.Meas.NSetMeas * prot.Meas.NRepMeas;

        mat_ro{k} = prot.Meas.NImageCols;
        mat_pe{k} = prot.Meas.NImageLins;
        dist_fact{k} = 0;
        
        b3D = or(strcmp(prot.MeasYaps.sKSpace.ucDimension,'0x04'),strcmp(prot.MeasYaps.sKSpace.ucDimension,'0x4'));
        if b3D
            mat_3D{k} = prot.MeasYaps.sKSpace.lPartitions;
            if isfield(prot.MeasYaps.sKSpace,'dSliceOversamplingForDialog')
                mat_3D{k} = mat_3D{k}/(1+prot.MeasYaps.sKSpace.dSliceOversamplingForDialog);
            end
        else
            mat_3D{k} = 1;
            if isfield(prot.MeasYaps.sGroupArray.asGroup{1},'dDistFact')
                dist_fact{k} = prot.MeasYaps.sGroupArray.asGroup{1}.dDistFact;
            end
        end
        mat_sl{k} = prot.Meas.NSlc;

        fov_ro{k} = prot.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
        fov_pe{k} = prot.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;
        fov_sl{k} = prot.MeasYaps.sSliceArray.asSlice{1}.dThickness;
        
        if b3D
            fov{k} = sprintf('%.fx%.fx%.f',fov_ro{k},fov_pe{k},fov_sl{k});
        else
            fov{k} = sprintf('%.fx%.fx%.1f',fov_ro{k},fov_pe{k},fov_sl{k});
        end

        slices{k} = prot.Meas.NSlc;
        
        res_ro{k} = fov_ro{k}/mat_ro{k};
        res_pe{k} = fov_pe{k}/mat_pe{k};
        res_sl{k} = fov_sl{k}/mat_3D{k};
        resolution{k} = sprintf('%.2fx%.2fx%.2f',res_ro{k},res_pe{k},res_sl{k});
        
        R_pe{k} = prot.MeasYaps.sPat.lAccelFactPE;
        R_3d{k} = prot.MeasYaps.sPat.lAccelFact3D;

        pf_pe{k} = evalPFmode(prot.MeasYaps.sKSpace.ucPhasePartialFourier);
        pf_3d{k} = evalPFmode(prot.MeasYaps.sKSpace.ucSlicePartialFourier);

        if isfield(prot.MeasYaps.sKSpace,'ucEnableEllipticalScanning')
            bEllipticScan{k} = strcmp(prot.MeasYaps.sKSpace.ucEnableEllipticalScanning,'0x1');
        else
            bEllipticScan{k} = 0;
        end
        
        notes{k} = '';
        if TI{k}~=0
            notes{k} = addnote(notes{k},sprintf('TI=%.3f s',TI{k}));
        end
        
        if numel(TR2{k})>0
            notes{k} = addnote(notes{k},sprintf('TR2=%.3f s',TR2{k}));
        end
        
        if ~b3D
            notes{k} = addnote(notes{k},sprintf('sl.dist=%.f%%',100*dist_fact{k}));
        end
        
        if R_pe{k}*R_3d{k}>1
            notes{k} = addnote(notes{k},strcat('R=',num2str(R_pe{k}),'x',num2str(R_3d{k})));
        end
        
        if ~strcmp(pf_pe{k},'off')
            notes{k} = addnote(notes{k},strcat(pf_pe{k},'y'));
        end
        
        if ~strcmp(pf_3d{k},'off')
            notes{k} = addnote(notes{k},strcat(pf_3d{k},'z'));
        end
        
        if bEllipticScan{k}
            notes{k} = addnote(notes{k},'ellipt.');
        end 
        
    end
    seq = seq.';
    mid = mid.';
    resolution = resolution.';
    fov = fov.';
    slices = slices.';
    voltage = voltage.';
    FA = FA.';
    TR = TR.';
    TE = TE.';
    TRvol = TRvol.';
    rep = rep.';
    notes = notes.';
    tab = table(seq,mid,resolution,fov,slices,FA,voltage,TR,TE,TRvol,rep,notes);
end

function notes = addnote(notes,note)
    if numel(notes)>0
        notes = [notes '; ' note];
    else
        notes = note;
    end
end

function out = evalPFmode(mode)
    switch mode
        case {'0x1','0x01'}
            out = '4/8';
        case {'0x2','0x02'}
            out = '5/8';
        case {'0x4','0x04'}
            out = '6/8';
        case {'0x8','0x08'};
            out = '7/8';
        case {'0x10'}
            out = 'off';
        case {'0x20'}
            out = 'auto';
        otherwise
            out = 'error';
    end
end
