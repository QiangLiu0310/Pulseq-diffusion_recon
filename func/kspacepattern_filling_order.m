function KspaceMaskKyKzt = kspacepattern_filling_order (dims,pat,varargin )
% dims are the dimensions of the of the image to be acquired x,y,z,t
% pat is a kernal sampling pattern in ky,kz
% the third argument is shift to be performed on the echo time to introduce some incoerence
% across echoes

% keyboard
dimspat = [1 size(pat)];

% initialize the shift to be performed along echo times
if nargin ==2
    shift = [1 1]
else
    shift = varargin{1};
end;

% ensure the imput dims have 4 dimensions
if length(size(dims))==3
    dims(end+1) = 1;
end

visualizeshift = 1;
if visualizeshift == 1       
figure();
end

KspaceMaskKyKzt = zeros(dims(2:4));

for nte = 1:dims(4)
    patmask = pat;
    patmask (pat~=1 )=0;
    patmask = circshift(patmask,shift * nte);
    patmaskorder = patmask;
%     patmaskorder (patmask==1) = 1:length(patmask==1);
    patmaskorder (patmask==1) = 1:numel(find(patmask==1));
    NpointsPerKernel = length(patmask==1);
    
    % KspaceMaskKyKzt = repmat(patmask,dims(2:3)./dimspat(2:3));
    if visualizeshift == 1       
        subplot(floor(sqrt(dims(4))),ceil(dims(4)/floor(sqrt(dims(4)))),nte)        
        imab(patmask) ; 
        title (['pattern of echo ',num2str(nte)] )
    end
    
    
    counter = 0;
    for kyblock = 1: dims(2)/dimspat(2)
        for kzblock = 1: dims(3)/dimspat(3)
            counter = counter + 1;
            KspaceMaskKyKzt((kyblock-1)*dimspat(2)+1:(kyblock)*dimspat(2),...
                (kzblock-1)*dimspat(3)+1:(kzblock)*dimspat(3),nte) = ...
                patmask.*(patmaskorder +(counter-1)*NpointsPerKernel) ;
        end;
    end;
end;

