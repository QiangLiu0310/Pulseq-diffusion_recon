function [K, w] = MultisliceGRAPPA_2kernal_leakBlock(K_Collapsed, varargin)

%K = MultisliceGRAPPA_2kernal...(K_Collapsed,K_Indiv,KernalSize,DataForFit,prot) :
%K = MultisliceGRAPPA_2kernal...(K_Collapsed,w,KernalSize)

% October 24 2010
% Kawin Setsompop

% add leak block (sp-SG) kernel calc to reduce leakage artifact
% Kawin Setsompop
% June 2012

% In this method, two Grappa kernels are created and applied, one for the odd line of k-space and one for the even line of k-space.
% This significantly reduce ghosting artifact. 
% As described in Moeller et al ghost correction can be an issue with multislice data as only the 'average ghost'
% can be corrected in the collapsed data. Moeller suggested a second step ghost correction (after the slice seperation) to mitigate this issue....this can help a bit but 
% the ghost issue really need to be accounted for during the grappa application process to really get good result. 
% The collapsed data k-space lines are not aligned due to imperfect ghost correction. By having a different grappa kernel for the odd and the even line this issue can be overcome. 
%
% Inputs:
% K_Collapsed  : the slice collapsed k space data 
% K_Indiv      : the 'acs' data. the data from individually acquired slice (for one simultanous multislice group)
% KernelSize   : [ NumberOfElementsAlongKx NumberOfElementsAlongKy distanceBetweenElementsAlongKx distanceBetweenElementsAlongdKy]
% DataForFit   : Part of the acs data that will be use in grappa kernel generation 
%                This can be specify 2 ways:
%                DataForFit = 'full'; use all the data from the acs dataset to estimate the grappa kernel
%                DataForFit = [SizeAlongKx SizeAlongKy]; use the center part of the acs k-space of size [SizeAlongKx SizeAlongKy] to estimate the kernel
%
% Outputs:
% K            : the slice seperated dataset
% w            : kernels 
% 
% K space data matrix dimensions are based on standard siemen convention
%   [#01]  NColMeas: 
%   [#02]  NLinMeas: 
%   [#03]  NChaMeas: 
%   [#04]  NSetMeas: 
%   [#05]  NEcoMeas: 
%   [#06]  NPhsMeas: 
%   [#07]  NRepMeas: 
%   [#08]  NSegMeas: 
%   [#09]  NParMeas: 
%   [#10]  NSlcMeas: 
%   [#11]  NIdaMeas: 
%   [#12]  NIdbMeas: 
%   [#13]  NIdcMeas: 
%   [#14]  NIddMeas: 
%   [#15]  NIdeMeas: 
%   [#16]  NAveMeas: 


debug = 0;

solver = 'pinv'; chi = []; eta = []; verbose = 0;

Nchannels = size(K_Collapsed,3);

if size(varargin,2) <= 3 % i.e. w is alredy calculated just apply
    w = varargin{1};
    KernalSize = varargin{2};
    KernalSupplied = 1;
    NslicesEX = size(w.odd,3);
    if size(varargin,2) == 3
        debug = varargin{3};
    end
else
    KernalSupplied = 0;
    K_Indiv = varargin{1};
    KernalSize = varargin{2};
    DataForFit = varargin{3};    
    prot = varargin{4};
    if size(varargin,2) == 5
        debug = varargin{5};
    end
    NslicesEX = max(size(K_Indiv));

    if isfield(prot,'sSolver')
      solver = prot.sSolver;
    end;
    if isfield(prot,'dChi')
      chi = prot.dChi;
    end;
    if isfield(prot,'dEta')
      eta = prot.dEta;
    end;
end

if length(KernalSize) == 2 % dim 3 and 4 are for SkipSteps
    KernalSize(3) = 1;
    KernalSize(4) = 1;
end

if KernalSupplied == 0
    s = size(K_Indiv{1});
%     K_CollapsedFake = zeros(s);
%     for count = 1:NslicesEX
%         K_CollapsedFake = K_CollapsedFake + K_Indiv{count};
%     end
    
    %% create Indexing
    if isfloat(DataForFit) % use only part of the supplied ACS data
        Nread = size(K_Collapsed,1);
        Nlin_full = prot.lPhaseEncodingLines;
        Nlin_part = size(K_Collapsed,2);
        
        DataEdgeBegin = [1 1] + floor(([Nread Nlin_full] - DataForFit)/2);
        DataEdgeEnd   = [Nread Nlin_full] - ceil(([Nread Nlin_full] - DataForFit)/2);
        
        if Nlin_full ~= Nlin_part % i.e. Partial Fourier so need to make sure use data in the center area (around dc) for kernal fitting
            pad_lines = Nlin_full - Nlin_part;
            DataEdgeBegin(2) = DataEdgeBegin(2) - pad_lines; % assuming missing part is 'pre'
            DataEdgeEnd(2) = DataEdgeEnd(2) - pad_lines;
            if DataEdgeBegin(2) < 1
                display('fitting data area exceed acquired data because of P.F')
                keyboard
            end
        end
    else
        % use full ACS data
        DataForFit = [size(K_Collapsed,1) size(K_Collapsed,2)];
        DataEdgeBegin = [1 1] ;
        DataEdgeEnd   = DataForFit;    
    end
    
    KxIndex = -floor( KernalSize(1)/2):floor((KernalSize(1)-1)/2); % i.e. kernal will be bulking in the "up-left" direction if kernal_size is even
    KyIndex = -floor( KernalSize(2)/2):floor((KernalSize(2)-1)/2);   
    KxIndex = KxIndex*KernalSize(3);
    KyIndex = KyIndex*KernalSize(4);
        
    BeginIndex = DataEdgeBegin  + floor(KernalSize(1:2)/2).*[KernalSize(3) KernalSize(4)]; %current mid point of the kernal
    EndIndex = DataEdgeEnd - floor((KernalSize(1:2)-1)/2).*[KernalSize(3) KernalSize(4)];
    
         
    %% Fitting Kernal
    
    % Kernal for ODD
    
    %Calculate DataMatrix
    tic
  
    MatrixSizeSliceSeperation = [(DataForFit(1) - (KernalSize(1)-1)*KernalSize(3))*ceil( (DataForFit(2) - (KernalSize(2)-1)*KernalSize(4)) /2 ) , ...     
                                    ( KernalSize(1) * KernalSize(2) )*Nchannels ];
    MatrixSizeSliceSeperation_leakBlock = [ MatrixSizeSliceSeperation(1)*NslicesEX, MatrixSizeSliceSeperation(2)];
    % tell you how well condition the problem is i.e. get good Kernal weights if m>>n where [m,n] = size(DataMatrix)
    
    for SliceCount = 1:NslicesEX
        DataMatrixCurrent = zeros(MatrixSizeSliceSeperation);
        
        CurrentIndex = BeginIndex;
        count = 1;
        
        for countY = 1:2:(DataForFit(2) - (KernalSize(2)-1)*KernalSize(4))
            for countX = 1:(DataForFit(1) - (KernalSize(1)-1)*KernalSize(3))
                DataMatrixCurrent(count,:) = reshape( K_Indiv{SliceCount}(KxIndex+CurrentIndex(1),KyIndex+CurrentIndex(2),:) ,1,[],1);
                count = count+1;
                CurrentIndex = CurrentIndex + [1 0];
            end
            CurrentIndex(1) = BeginIndex(1);
            CurrentIndex(2) = CurrentIndex(2) + 2;
        end

        if ( ~isempty(eta) && (eta ~= 0) )
          % compute the power in each line of the linsys
          p_src{SliceCount} = sum(abs(DataMatrixCurrent).^2,2);

          % and apply the inverse weighting to _normalize_ the linear system
          p_src{SliceCount}(find( p_src{SliceCount} )) = p_src{SliceCount}(find( p_src{SliceCount} )).^(-eta/2);
          DataMatrixCurrent = diag(p_src{SliceCount}) * DataMatrixCurrent ;
        end;
        
        if SliceCount == 1
            DataMatrix = DataMatrixCurrent;
        else
            DataMatrix = cat(1,DataMatrix,DataMatrixCurrent);
        end
        
    end
    
    disp_debug_msg(['ODD: DataMatrix Form: ' num2str(toc) ' s'],debug);
    
    % Calculate Weights
    tic
    %cond(DataMatrix,2)

    if strcmp(solver,'pinv');
      InvDataMatrix = pinv(DataMatrix);
      disp_debug_msg(['Inversion: ' num2str(toc) ' s'],debug);
    else
      AA = DataMatrix'*DataMatrix;
      At = DataMatrix';
    end;

    flag = 1;
    % figure; imagesc(abs(DataMatrix)); 
    for ChCount = 1:Nchannels
        for SliceCount = 1:NslicesEX
            ACS_actual = K_Indiv{SliceCount}(BeginIndex(1):EndIndex(1),BeginIndex(2):2:EndIndex(2),ChCount);
            if flag == 1
                sV = size(ACS_actual(:));
                flag = 0;
            end
            if ( ~isempty(eta) && (eta~=0) )
              ACS_actual(:) = ACS_actual(:) .* p_src{SliceCount};
            end;
            ACS = cat(1,zeros(sV(1)*(SliceCount-1),1),ACS_actual(:), zeros(sV(1)*(NslicesEX -SliceCount),1));
            if strcmp(solver,'pinv');
              w.odd(:,ChCount,SliceCount) = InvDataMatrix*ACS(:);
            elseif strcmp(solver,'cgsr');
              w.odd(:,ChCount,SliceCount) = cgsolv( AA + chi*eye(size(AA)), At*ACS(:), zeros(size(At,1),1), size(At,1) );
            else
              w.odd(:,ChCount,SliceCount) = bicg( AA + chi*eye(size(AA)), At*ACS(:) );
            end;
        end
    end
   
    % Kernal for EVEN
    
    %Calculate DataMatrix
    tic
    
    MatrixSizeSliceSeperation = [(DataForFit(1) - (KernalSize(1)-1)*KernalSize(3))*floor( (DataForFit(2) - (KernalSize(2)-1)*KernalSize(4)) /2 ) , ...
        ( KernalSize(1) * KernalSize(2) )*Nchannels ];
    MatrixSizeSliceSeperation_leakBlock = [ MatrixSizeSliceSeperation(1)*NslicesEX, MatrixSizeSliceSeperation(2)];
    % tell you how well condition the problem is i.e. get good Kernal weights if m>>n where [m,n] = size(DataMatrix)
    
    for SliceCount = 1:NslicesEX
        DataMatrixCurrent = zeros(MatrixSizeSliceSeperation);
        
        CurrentIndex = BeginIndex;
        CurrentIndex(2) = CurrentIndex(2)+1;
        count = 1;
        for countY = 2:2:(DataForFit(2) - (KernalSize(2)-1)*KernalSize(4))
            for countX = 1:(DataForFit(1) - (KernalSize(1)-1)*KernalSize(3))
                DataMatrixCurrent(count,:) = reshape(K_Indiv{SliceCount}(KxIndex+CurrentIndex(1),KyIndex+CurrentIndex(2),:) ,1,[],1);
                count = count+1;
                CurrentIndex = CurrentIndex + [1 0];
            end
            CurrentIndex(1) = BeginIndex(1);
            CurrentIndex(2) = CurrentIndex(2) + 2;
        end
        if ( ~isempty(eta) && (eta ~= 0) )
          % compute the power in each line of the linsys
          p_src{SliceCount} = sum(abs(DataMatrixCurrent).^2,2);

          % and apply the inverse weighting to _normalize_ the linear system
          p_src{SliceCount}(find( p_src{SliceCount} )) = p_src{SliceCount}(find( p_src{SliceCount} )).^(-eta/2);
          DataMatrixCurrent = diag(p_src{SliceCount}) * DataMatrixCurrent ;
        end;
        
        if SliceCount == 1
            DataMatrix = DataMatrixCurrent;
        else
            DataMatrix = cat(1,DataMatrix,DataMatrixCurrent);
        end
        
    end
        
        
    disp_debug_msg(['EVEN: DataMatrix Form: ' num2str(toc) ' s'],debug);
        
    % Calculate Weights
    tic
    if strcmp(solver,'pinv');
      InvDataMatrix = pinv(DataMatrix);
      disp_debug_msg(['Inversion: ' num2str(toc) ' s'],debug);
    else
      AA = DataMatrix'*DataMatrix;
      At = DataMatrix';
      if verbose, keyboard; end;
    end;

    disp_debug_msg(['Inversion: ' num2str(toc) ' s'],debug);
    
    flag = 1;
    % figure; imagesc(abs(DataMatrix));
    for ChCount = 1:Nchannels
        for SliceCount = 1:NslicesEX
            ACS_actual = K_Indiv{SliceCount}(BeginIndex(1):EndIndex(1),BeginIndex(2)+1:2:EndIndex(2),ChCount);
            if flag == 1
                sV = size(ACS_actual(:));
                flag = 0;
            end
            if ( ~isempty(eta) && (eta~=0) )
              ACS_actual(:) = ACS_actual(:) .* p_src{SliceCount};
            end;
            ACS = cat(1,zeros(sV(1)*(SliceCount-1),1),ACS_actual(:), zeros(sV(1)*(NslicesEX -SliceCount),1));
            if strcmp(solver,'pinv');
              w.even(:,ChCount,SliceCount) = InvDataMatrix*ACS(:);
            elseif strcmp(solver,'cgsr');
              w.even(:,ChCount,SliceCount) = cgsolv( AA + chi*eye(size(AA)), At*ACS(:), zeros(size(At,1),1), size(At,1) );
            else
              w.even(:,ChCount,SliceCount) = bicg( AA + chi*eye(size(AA)), At*ACS(:) );
            end;
        end
    end

else
    % need to determine which set of grappa kernels is for odd and which is
    % for even lines.
    
    % assume using the whole dataset during the fitting.....
    
    DataForFit = [size(K_Collapsed,1) size(K_Collapsed,2)];
    DataEdgeBegin = [1 1] ;
    DataEdgeEnd   = DataForFit;
    
    BeginIndex = DataEdgeBegin  + floor(KernalSize(1:2)/2).*[KernalSize(3) KernalSize(4)]; %current mid point of the kernal
    EndIndex = DataEdgeEnd - floor((KernalSize(1:2)-1)/2).*[KernalSize(3) KernalSize(4)];
end


%% Seperate Slices using the calculated Kernal

METHOD = 'conv2';
tic

switch METHOD,
    %----------------------------------------------------------------------%
    case 'conv2',
        s = size(K_Collapsed);
       if length(s) < 10
           s(end+1:9) = 1;
       end
       
       s(10) = NslicesEX;
       K_odd = zeros(s);
       K_even = zeros(s);
       
       Wcurrent_odd = zeros((KernalSize(1)-1)*KernalSize(3)+1,(KernalSize(2)-1)*KernalSize(4)+1,Nchannels);
       Wcurrent_even = zeros((KernalSize(1)-1)*KernalSize(3)+1,(KernalSize(2)-1)*KernalSize(4)+1,Nchannels);
       
       for SliceCount = 1:NslicesEX
           for ChCount = 1:Nchannels
               Wcurrent_odd(1:KernalSize(3):end,1:KernalSize(4):end,:)  = reshape(w.odd(:,ChCount,SliceCount),KernalSize(1),KernalSize(2),Nchannels);
               Wcurrent_even(1:KernalSize(3):end,1:KernalSize(4):end,:)  = reshape(w.even(:,ChCount,SliceCount),KernalSize(1),KernalSize(2),Nchannels);
               for SetCount = 1:size(K_Collapsed,4)
                   for AvgCount = 1:size(K_Collapsed,7)
                       for ChCount2 = 1:Nchannels
                           WcurrentCh2_odd = Wcurrent_odd(end:-1:1,end:-1:1,ChCount2);
                           WcurrentCh2_even = Wcurrent_even(end:-1:1,end:-1:1,ChCount2);
                           K_odd(:,:,ChCount,SetCount,:,:,AvgCount,:,:,SliceCount) =  K_odd(:,:,ChCount,SetCount,:,:,AvgCount,:,:,SliceCount) + conv2(K_Collapsed(:,:,ChCount2,SetCount,:,:,AvgCount,:,:,:),WcurrentCh2_odd,'same') ;
                           K_even(:,:,ChCount,SetCount,:,:,AvgCount,:,:,SliceCount) =  K_even(:,:,ChCount,SetCount,:,:,AvgCount,:,:,SliceCount) + conv2(K_Collapsed(:,:,ChCount2,SetCount,:,:,AvgCount,:,:,:),WcurrentCh2_even,'same') ;
                       end
                   end
               end
           end
       end
       
       if rem(BeginIndex(2),2) == 1
           disp_debug_msg('first Kernal fit odd lines well',debug);
           disp_debug_msg(' ',debug);
           K = zeros(s);
           K(:,1:2:end,:,:,:,:,:,:,:,:) = K_odd(:,1:2:end,:,:,:,:,:,:,:,:);
           K(:,2:2:end,:,:,:,:,:,:,:,:) = K_even(:,2:2:end,:,:,:,:,:,:,:,:);
       else
           disp_debug_msg('first Kernal fit even lines well',debug);
           disp_debug_msg(' ',debug);
           K = zeros(s);
           K(:,1:2:end,:,:,:,:,:,:,:,:) = K_even(:,1:2:end,:,:,:,:,:,:,:,:);
           K(:,2:2:end,:,:,:,:,:,:,:,:) = K_odd(:,2:2:end,:,:,:,:,:,:,:,:);
       end
    case 'ImageDomain'
        
        PadKx = round(KernalSize(1)*KernalSize(3)/2);
        PadKy = round(KernalSize(2)*KernalSize(4)/2);
        
        I_Collapsed_Padded = mrir_iDFT(mrir_iDFT(padarray(K_Collapsed,[PadKx PadKy]),1),2); 
        dims = size(I_Collapsed_Padded);
        if length(dims) < 10
            dims(end+1:9) = 1;
        end
        dims(10) = NslicesEX;
        
        K_odd = zeros(dims);
        K_even = zeros(dims);
        K_Kernel_odd = zeros([dims(1:3), ones(1,6), dims(10)]);
        K_Kernel_even = zeros([dims(1:3), ones(1,6), dims(10)]);
        
        KxIndex = -floor( KernalSize(1)/2):floor((KernalSize(1)-1)/2); % i.e. kernal will be bulking in the "up-left" direction if kernal_size is even
        KyIndex = -floor( KernalSize(2)/2):floor((KernalSize(2)-1)/2);
        %KxIndex = KxIndex*KernalSize(3) + ceil((dims(1)+1)/2);
        %KyIndex = KyIndex*KernalSize(4) + ceil((dims(2)+1)/2);
        KxIndex = KxIndex*KernalSize(3) + floor((dims(1)+1)/2); % this is because convolution definition require the flipping of the kernel (and move center of k-spce for fft) so instead of ceil here need to do floor
        KyIndex = KyIndex*KernalSize(4) + floor((dims(2)+1)/2);
        
        
        for ChCount = 1:Nchannels
            K_Kernel_odd(KxIndex,KyIndex,:,:,:,:,:,:,:,:)  = reshape(w.odd(:,ChCount,:),[KernalSize(1),KernalSize(2),Nchannels,1,1,1,1,1,1,NslicesEX]);
            I_Kernel_odd = mrir_iDFT(mrir_iDFT(K_Kernel_odd(end:-1:1,end:-1:1,:,1,1,1,1,1,1,:),1),2);
            K_Kernel_even(KxIndex,KyIndex,:,:,:,:,:,:,:,:)  = reshape(w.even(:,ChCount,:),[KernalSize(1),KernalSize(2),Nchannels,1,1,1,1,1,1,NslicesEX]);
            I_Kernel_even = mrir_iDFT(mrir_iDFT(K_Kernel_even(end:-1:1,end:-1:1,:,1,1,1,1,1,1,:),1),2);
            
            K_odd(:,:,ChCount,:,:,:,:,:,:,:) = sum(mrir_fDFT(mrir_fDFT(repmat(I_Kernel_odd,[1 1 1 dims(4) 1 1 dims(7) 1 1 1]).*repmat(I_Collapsed_Padded,[1 1 1 1 1 1 1 1 1 NslicesEX]) ,1),2),3);
            K_even(:,:,ChCount,:,:,:,:,:,:,:) = sum(mrir_fDFT(mrir_fDFT(repmat(I_Kernel_even,[1 1 1 dims(4) 1 1 dims(7) 1 1 1]).*repmat(I_Collapsed_Padded,[1 1 1 1 1 1 1 1 1 NslicesEX]) ,1),2),3);
%             for SetCount = 1:size(K_Collapsed,4)
%                 for AvgCount = 1:size(K_Collapsed,7)
%                     K_odd(:,:,ChCount,SetCount,:,:,AvgCount,:,:,:) = sum(mrir_fDFT(mrir_fDFT(I_Kernel_odd.*repmat(I_Collapsed_Padded(:,:,:,SetCount,:,:,AvgCount,:,:,:),[1 1 1 1 1 1 1 1 1 NslicesEX]) ,1),2),3);
%                     K_even(:,:,ChCount,SetCount,:,:,AvgCount,:,:,:) = sum(mrir_fDFT(mrir_fDFT(I_Kernel_even.*repmat(I_Collapsed_Padded(:,:,:,SetCount,:,:,AvgCount,:,:,:),[1 1 1 1 1 1 1 1 1 NslicesEX]) ,1),2),3);
%                 end
%             end
        end
        
        dims(1) = size(K_Collapsed,1);
        dims(2) = size(K_Collapsed,2);
        if rem(BeginIndex(2),2) == 1
            disp('first Kernal fit odd lines well')
            disp(' ')
            K = zeros(dims);
            K(:,1:2:end,:,:,:,:,:,:,:,:) = K_odd(1+PadKx:end-PadKx,1+PadKy:2:end-PadKy,:,:,:,:,:,:,:,:);
            K(:,2:2:end,:,:,:,:,:,:,:,:) = K_even(1+PadKx:end-PadKx,2+PadKy:2:end-PadKy,:,:,:,:,:,:,:,:);
        else
            disp('first Kernal fit even lines well')
            disp(' ')
            K = zeros(dims);
            K(:,1:2:end,:,:,:,:,:,:,:,:) = K_even(1+PadKx:end-PadKx,1+PadKy:2:end-PadKy,:,:,:,:,:,:,:,:);
            K(:,2:2:end,:,:,:,:,:,:,:,:) = K_odd(1+PadKx:end-PadKx,2+PadKy:2:end-PadKy,:,:,:,:,:,:,:,:);
        end
end

disp_debug_msg(['Kernal Application: ' num2str(toc) ' s'],debug);
%% display

if debug == 1
    Set = 1; Avg = 1;
    I = mrir_array_combine(mrir_conventional_2d(K(:,:,:,Set,:,:,Avg,:,:,:)),0);

    figure(1); imagesc(mrir_array_combine(mrir_conventional_2d(K_Collapsed(:,:,:,Set,:,:,Avg,:,:,:)),0)); 
    axis equal; title('Acquired Collasped Image'); colormap gray;
    %figure(2); subplot(2,NslicesEX,1); imagesc(abs(K_Collapsed(:,:,1))); title('KspaceCollapsedForChn1') % plot result of kspace fitting for channel 1
    
    for SliceCount = 1:NslicesEX
        figure(3);  
        subplot(2,NslicesEX,SliceCount); imagesc(mrir_array_combine(mrir_conventional_2d(K_Indiv{SliceCount}(:,:,:,Set,:,:,Avg,:,:,:)),0));
        axis equal; title(['Indiv Slice' num2str(SliceCount) ]); colorbar; colormap gray;
        subplot(2,NslicesEX,NslicesEX+SliceCount); imagesc(I(:,:,SliceCount));
        axis equal; title(['BackOut Slice' num2str(SliceCount) ]); colorbar; colormap gray;
    end
    keyboard
end


function disp_debug_msg( str, flag )

if (flag)
  disp(str);
end;