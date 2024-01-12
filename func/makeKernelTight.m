function [pat ] = makeKernelTight(Rx, Ry, delta)
% function [pat stepSize startPos] = makeKernelJose(Rx, Ry, delta,varargin)
% creates a grappa kernel with acceleration Rx x Ry and a CAIPI step of
% delta
% it also outputs the StepSize and startPos needed for the grappa code to
% populate the whole k-space

visualize=0;

%% defining the pat matrix

pat = zeros( Ry * (Rx), Rx  * (Ry) ); % creates the Rx * Ry
% source points are marked with 1
% target points are marked with 0.5
pat(1:Rx,1:Ry ) = 0.5;
pat(1:Rx:end,1:Ry:end) = 1;


%% caipi kernel with aproximately the correct number of source points
patStarts = zeros(Ry*Rx, Rx*Ry );
patStarts( 1:Rx:end,1:Ry:end) = 1;

%% visualizing  comment out
if visualize == 1
    figure(1)
    subplot(221)
    imab(pat(:,:,1))
end

%% introduces the caipi in the small 2x2 kernel
dims= size(pat);
l=0;
for k=1:Rx:dims(1)
    l=l+1;
    pat((l-1)*Rx+[1:Rx],:)=circshift(pat((l-1)*Rx+[1:Rx],:),[0 (l-1)*delta]);
end;

if Rx*Ry>1
    startPosind=find(pat==0.5);
    [coord(1,:) coord(2,:) ] =ind2sub(dims,startPosind);
    
    %% introduces the caipi delta in the Rx*Ry x Rx*Ry which allows computing the start positions
    dims= size(patStarts);
    l=0;
    for k=1:Rx:dims(1)
        l=l+1;
        patStarts((l-1)*Rx+[1:Rx],:)=circshift(patStarts((l-1)*Rx+[1:Rx],:),[0 (l-1)*delta]);
    end;
    stepSize=dims;
    startPosind=find(patStarts);
    [startPos(1,:) startPos(2,:) ] =ind2sub(dims,startPosind);
    
    %% search of the best targets within the kernel
    
    % positions has the information of the coordinates of the targets
    % associeted with each source point in the Rx,Ry matrix
    pos=bsxfun(@plus,startPos-1,reshape(coord,[2 1 Rx*Ry-1]));
    pos(2,:,:)=mod(pos(2,:,:)-1,dims(2))+1;
    
    % calculates the distance to of the eache target to all the sources and
    % assumes the correct target is the one that has the minimum distance
    distToAllSources = sum(sum(bsxfun(@minus,pos,reshape(startPos,[2 1 1 size(startPos,2)])).^2,4),1);
    [temp temppos] = min(distToAllSources,[],2);
    
    % populates the targets
    for k=1:length(temppos)
        patStarts(pos(1,temppos(k),k),pos(2,temppos(k),k))=0.5;
    end;
    
    pat=patStarts;
end
if visualize==1
    subplot(223)
    imab(patStarts)
    sources=length(find(patStarts==0.5));
    targets=length(find(patStarts==1));
    EffectiveKernel=sqrt((targets+1)*sources/(Rx*Ry));
    ylabel(['Rx',num2str(Rx),'Ry',num2str(Ry),'\Delta',num2str(delta),'Kernel~',num2str(round(EffectiveKernel))])
end;

% keyboard
