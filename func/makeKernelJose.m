function [pat stepSize startPos] = makeKernelJose(Rx, Ry, delta,varargin)
% function [pat stepSize startPos] = makeKernelJose(Rx, Ry, delta,varargin)
% creates a grappa kernel with acceleration Rx x Ry and a CAIPI step of
% delta
% it also outputs the StepSize and startPos needed for the grappa code to
% populate the whole k-space

visualize=0;

%% defining the pat matrix

% create a standard grappa kernel
if nargin==3;
    nx = 3; ny = 3;
    nx_small = 2; ny_small = 2;
elseif nargin==4
    nx = varargin{1}; ny = varargin{1};
    nx_small = 2; ny_small = 2;
end;

% small 2x2 kernel
pat = zeros( 2 * (Rx), 2  * (Ry) );
start_x = min(2 + (floor(nx_small/2) - 1) * Rx,Rx*nx_small);
start_y = min(2 + (floor(ny_small/2) - 1) * Ry,Ry*ny_small);

% source points are marked with 1
% target points are marked with 0.5
pat(start_x:start_x+Rx-1,start_y:start_y+Ry-1 ) = 0.5;
pat( 1:Rx:end,1:Ry:end) = 1;

%% the default size of the the caipi matrix is (Rx*Ry)^2 and it will be first shrinked to the minimum kernel
%  and then grown until it is aproximately the nxn kernel size

if mod (Ry,delta) == 0;
    Ry_new = Ry/delta;
else
    if mod (Ry,delta) == Ry;
        Ry_new = 1;
    else
        Ry_new = Ry;
    end
end

%% caipi kernel with aproximately the correct number of source points
% patStarts = zeros(ceil(nx/Ry_new)*Ry_new*Rx, ceil(ny/Rx)*Rx*Ry );
patStarts = zeros(ceil(nx/Ry_new)*Ry_new*Rx, ceil(ny)*Ry );
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
startPosind=find(pat==0.5);
[coord(1,:) coord(2,:) ] =ind2sub(dims,startPosind);

%% visualizing, if possible, the version from mayur for comparison
if visualize==1    
    subplot(222)
    imab(pat(:,:,1))    
    try
        [kernel stepSizeb startPosb] = makeKernel(Rx, Ry, delta);
        subplot(224)        
        imab(kernel),
        sources=length(find(kernel==0.5));
        targets=length(find(kernel==1));
        EffectiveKernel=sqrt((targets+1)*sources/(Rx*Ry));
        ylabel(['Rx',num2str(Rx),'Ry',num2str(Ry),'\Delta',num2str(delta),'Kernel~',num2str(round(EffectiveKernel))])        
        title('Hard coded Mayur');
        stepSizeb;
        startPosb;
    end;
end


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

if visualize==1    
    subplot(223)
    imab(patStarts)
    sources=length(find(patStarts==0.5));
    targets=length(find(patStarts==1));
    EffectiveKernel=sqrt((targets+1)*sources/(Rx*Ry));    
    ylabel(['Rx',num2str(Rx),'Ry',num2str(Ry),'\Delta',num2str(delta),'Kernel~',num2str(round(EffectiveKernel))])
end;

% keyboard
