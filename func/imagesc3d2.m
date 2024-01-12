function [  ] = imagesc3d2( img, pos, fig_num, rot_deg, scale_fig, avg_size, title_fig, MIP, axis_square )
%IMAGESC3D Summary of this function goes here
%   Detailed explanation goes here


img = single(img);

if norm(imag(img(:))) > 0
    img = abs(img);
end

if nargin < 2
    pos = round(size(img)/2);
else
    pos = round(pos);
end

if nargin < 3
    fig_num = 1;
end

if nargin < 4
    rot_deg = [0,0,0];
end

if nargin < 5
    scale_fig = [min(img(:)), max(img(:))];
end

if nargin < 6
    avg_size = 0;
end

if nargin < 7
    title_fig = '';
end

if nargin < 8
    MIP = 0;
end

if nargin < 9
    axis_square = 0;
end

Lavg_size = 0;
Ravg_size = 0;

if avg_size > 0
    if mod(avg_size,2)
        % odd
        Lavg_size = (avg_size-1)/2;
        Ravg_size = (avg_size-1)/2;
    else
        % even
        Lavg_size = (avg_size / 2) - 1;
        Ravg_size = (avg_size / 2);
    end
end


if MIP == 0
    if axis_square
        figure(fig_num), subplot(131), imagesc( imrotate(squeeze(mean(img(pos(1)-Lavg_size:pos(1)+Ravg_size,:,:),1)), rot_deg(1)), scale_fig ), axis square off, colormap gray
        figure(fig_num), subplot(132), imagesc( imrotate(squeeze(mean(img(:,pos(2)-Lavg_size:pos(2)+Ravg_size,:),2)), rot_deg(2)), scale_fig ), axis square off, colormap gray
            title(title_fig, 'color', 'g', 'fontsize', 24)
        figure(fig_num), subplot(133), imagesc( imrotate(mean(img(:,:,pos(3)-Lavg_size:pos(3)+Ravg_size),3), rot_deg(3)), scale_fig ), axis square off, colormap gray
    else
        figure(fig_num), subplot(131), imagesc( imrotate(squeeze(mean(img(pos(1)-Lavg_size:pos(1)+Ravg_size,:,:),1)), rot_deg(1)), scale_fig ), axis image off, colormap gray
        figure(fig_num), subplot(132), imagesc( imrotate(squeeze(mean(img(:,pos(2)-Lavg_size:pos(2)+Ravg_size,:),2)), rot_deg(2)), scale_fig ), axis image off, colormap gray
            title(title_fig, 'color', 'g', 'fontsize', 24)
        figure(fig_num), subplot(133), imagesc( imrotate(mean(img(:,:,pos(3)-Lavg_size:pos(3)+Ravg_size),3), rot_deg(3)), scale_fig ), axis image off, colormap gray
    end
else
    if MIP == 1
        if axis_square
            figure(fig_num), subplot(131), imagesc( imrotate(squeeze(min(img(pos(1)-Lavg_size:pos(1)+Ravg_size,:,:),[],1)), rot_deg(1)), scale_fig ), axis square off, colormap gray
            figure(fig_num), subplot(132), imagesc( imrotate(squeeze(min(img(:,pos(2)-Lavg_size:pos(2)+Ravg_size,:),[],2)), rot_deg(2)), scale_fig ), axis square off, colormap gray
                title(title_fig, 'color', 'g', 'fontsize', 24)
            figure(fig_num), subplot(133), imagesc( imrotate(min(img(:,:,pos(3)-Lavg_size:pos(3)+Ravg_size),[],3), rot_deg(3)), scale_fig ), axis square off, colormap gray    
        else
            figure(fig_num), subplot(131), imagesc( imrotate(squeeze(min(img(pos(1)-Lavg_size:pos(1)+Ravg_size,:,:),[],1)), rot_deg(1)), scale_fig ), axis image off, colormap gray
            figure(fig_num), subplot(132), imagesc( imrotate(squeeze(min(img(:,pos(2)-Lavg_size:pos(2)+Ravg_size,:),[],2)), rot_deg(2)), scale_fig ), axis image off, colormap gray
                title(title_fig, 'color', 'g', 'fontsize', 24)
            figure(fig_num), subplot(133), imagesc( imrotate(min(img(:,:,pos(3)-Lavg_size:pos(3)+Ravg_size),[],3), rot_deg(3)), scale_fig ), axis image off, colormap gray    
        end
    else
        if axis_square
            figure(fig_num), subplot(131), imagesc( imrotate(squeeze(max(img(pos(1)-Lavg_size:pos(1)+Ravg_size,:,:),[],1)), rot_deg(1)), scale_fig ), axis square off, colormap gray
            figure(fig_num), subplot(132), imagesc( imrotate(squeeze(max(img(:,pos(2)-Lavg_size:pos(2)+Ravg_size,:),[],2)), rot_deg(2)), scale_fig ), axis square off, colormap gray
                title(title_fig, 'color', 'g', 'fontsize', 24)
            figure(fig_num), subplot(133), imagesc( imrotate(max(img(:,:,pos(3)-Lavg_size:pos(3)+Ravg_size),[],3), rot_deg(3)), scale_fig ), axis square off, colormap gray    
        else
            figure(fig_num), subplot(131), imagesc( imrotate(squeeze(max(img(pos(1)-Lavg_size:pos(1)+Ravg_size,:,:),[],1)), rot_deg(1)), scale_fig ), axis image off, colormap gray
            figure(fig_num), subplot(132), imagesc( imrotate(squeeze(max(img(:,pos(2)-Lavg_size:pos(2)+Ravg_size,:),[],2)), rot_deg(2)), scale_fig ), axis image off, colormap gray
                title(title_fig, 'color', 'g', 'fontsize', 24)
            figure(fig_num), subplot(133), imagesc( imrotate(max(img(:,:,pos(3)-Lavg_size:pos(3)+Ravg_size),[],3), rot_deg(3)), scale_fig ), axis image off, colormap gray    
        end
    end
    
end

set(gcf, 'color', 'k')
set(gca, 'LooseInset', get(gca,'TightInset'))
    
end

