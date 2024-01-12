function Display_MainAtributes(seq, SeqParams);

[kfa,ta,kf,t,texcite]=seq.calculateKspacePP(); % ta is the timing of measurements and kfa the k-space coordinate at those times

TReff = mean(diff(texcite));
index = find(and((t>TReff*(SeqParams.Ndummy)),(t<TReff*(SeqParams.Ndummy+1))));
index2 = find(and((ta>TReff*(SeqParams.Ndummy)),(ta<TReff*(SeqParams.Ndummy+1))));

figure()
set(gcf,'Position', 1.0e+03 * [ -1.9177   -0.1703    1.4873    0.8113])

subplot(222)
time=(1:length(SeqParams.labelData0_Nav1))*SeqParams.TR;
hold off
plot( time , SeqParams.labelData0_Nav1,'.');
hold on
if ~isempty(find(SeqParams.labelData0_Nav1))
    timenav = linspace(SeqParams.TR,length(SeqParams.labelData0_Nav1)*SeqParams.TR,length(SeqParams.navcount));
    plot(timenav , mod(SeqParams.navcount,prod(SeqParams.RyRzCaipi(1:2)))*0.1+1,'r.')
    legend(['data (0) and nav (1) data'],['independent MC navs']);
    title([num2str(max(SeqParams.navcount)),'mc navigators acquires at 8mm res'])
    text(0,0.5, ['Fraction of samples in MC navigator spent fully sampling k-space \ centre ',num2str(length(find(SeqParams.labelData0_Nav1==1))/length(SeqParams.labelData0_Nav1))])
    text(0,0.2, ['currently each 8mm navigator acquired with acceleration factor ', num2str(prod(SeqParams.RyRzCaipi(1:2)))])
else
    legend(['data (0) and nav (1) data']);
    
end
xlabel(['duration of scan (s)'])
ylabel('Data vs Centre k-space')


subplot(231)
if isfield(SeqParams,'KspaceOrderNavKyKzt')
    measures = sum(double(SeqParams.KspaceMaskKyKzt~=0),3)+sum(double(SeqParams.KspaceOrderNavKyKzt~=0),3);
else
    measures = sum(double(SeqParams.KspaceMaskKyKzt~=0),3);
    
end
imab(measures),colorbar,colormap(gray);
xlabel('K_y'),ylabel('K_z'),axis square
title(['ntimes kykz is sampled over the ',num2str(length(SeqParams.TE)) ,' echo times'])

subplot(234)
imab(double(SeqParams.KspaceMaskKyKzt(:,:,1))),colorbar,colormap(gray);
xlabel('K_y'),ylabel('K_z'),axis square
title(['Shot number of used for kykz TE1'])
if length(SeqParams.TE)>1
    subplot(235)
    imab(double(SeqParams.KspaceMaskKyKzt(:,:,2))),colorbar,colormap(gray);
    xlabel('K_y'),ylabel('K_z'),axis square
    title(['Shot number of used for kykz TE2'])
end
% set(gcf,'Position',[2.6697    0.1043    0.5600    0.4200]*1000)

subplot(236)
plot3(kf(1,index),kf(2,index),kf(3,index));
hold on
plot3(kfa(1,index2),kfa(2,index2),kfa(3,index2),'.r');
axis tight
title ([ 'showing the k space trajectory after first excitation '] )

