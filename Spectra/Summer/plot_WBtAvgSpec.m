%% plot_WBtAvgSpec.m


    % Noah Clark
    % Created: 6/16/2023
    
    % Purpose: 
    %          - compute the average depths of each buoys
    %          - plot the time averaged spectrum for all wave buoys  
 
    
    
%%

clear;clc;

    %Load in data from WBvariables.mat:
load('WBvariables.mat')

%%

    %Determine the average depths of each buoy
for i = 1:3
    BavgD{i} = mean(Bdepth{i});
    XavgD{i} = mean(Xdepth{i});
end

    %China Rock (B):
figure(1);clf;
set(gcf,'position',[0,40,450,450])
for i = 1:3
    Bmeanspec{i} = nanmean(BSee{i}');  
    plot(Bfreq{i},Bmeanspec{i},plotcolors(i),'LineWidth',2)
    hold on
end
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Energy (m^2/Hz)')
    title('Time-Averaged Spectrum for the China Rock Buoys')
    xlim([0,0.5])
    ylim([0,1.2])
    legend(sprintf('B01 - Avg Depth (m): %g\n',round(BavgD{1},2)),...
        sprintf('B03 - Avg Depth (m): %g\n',round(BavgD{2},2)),...
        sprintf('B05 - Avg Depth (m): %g\n',round(BavgD{3},2)))

    %Asilomar (X):
figure(2);clf;
set(gcf,'position',[450,40,450,450])
for i = 1:3
    Xmeanspec{i} = nanmean(XSee{i}');
    plot(Xfreq{i},Xmeanspec{i},plotcolors(i),'LineWidth',2)
    hold on
end
    grid on
    xlabel('Frequency (Hz)')
    ylabel('Energy (m^2/Hz)')
    title('Time-Averaged Spectrum for the Asilomar Buoys')
    xlim([0,0.5])
    ylim([0,1.2])
    legend(sprintf('X01 - Avg Depth (m): %g\n',round(XavgD{1},2)),...
        sprintf('X03 - Avg Depth (m): %g\n',round(XavgD{2},2)),...
        sprintf('X04 - Avg Depth (m): %g\n',round(XavgD{3},2)))



%%

    % Saving Variables:
save('WBvariables.mat','-append','XavgD','BavgD','Xmeanspec','Bmeanspec')

