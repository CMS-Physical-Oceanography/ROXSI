%% plot_WBdirecto.m


    % Noah Clark
    % Created: 6/16/2023
    
    % Purpose: 
    %          - plot the EMEM spectrums("directograms") (wave direction 
    %            organized based on wave frequency and time) for each wave
    %            buoy   
    
    
%%

clear;clc;

    %Load in data from WBvariables.mat:
load('WBvariables.mat')


%%

LB1x1 = [datetime('02-Jul-2022 00:00:00','TimeZone','America/Los_Angeles')...
    , datetime('02-Jul-2022 00:00:00','TimeZone','America/Los_Angeles')];
LB1x2 = [datetime('06-Jul-2022 12:00:00','TimeZone','America/Los_Angeles')...
    , datetime('06-Jul-2022 12:00:00','TimeZone','America/Los_Angeles')];
LB1y1 = [0 0.55];
LB1y2 = [0 0.55];

%China Rock Buoys:
figure(1);clf;
set(gcf,'position',[800,300,500,700])
for i = 1:3
    subplot(3,1,i)
    pcolor(Btime{i},Gfreq,BEMEM{i})
    shading flat
    cb = colorbar;
    caxis([200 360])
    ylim([0 0.55])
    if i == 1
        hold on
        plot(LB1x1,LB1y1,'-r','LineWidth',2)
        plot(LB1x2,LB1y2,'-r','LineWidth',2)
        title({'China Rock Normal Wave Direction - 292.8^o N','','','B01 Directogram'})
    elseif i == 2
        title('B03 Directogram')
    else
        title('B05 Directogram')
        xlabel('Date/Time')
    end
    ylabel(cb,'Wave Direction (degrees)')
    ylabel('Frequency (Hz)')
end

%Asilomar Buoys:
figure(2);clf;
set(gcf,'position',[950,300,500,700])
for i = 1:3
    subplot(3,1,i)
    pcolor(Xtime{i},Gfreq,XEMEM{i})
    shading flat
    cb = colorbar;
    caxis([200 360])
    ylim([0 0.55])
    if i == 1
        title({'Asilomar Normal Wave Direction - 293.9^o N','','','X01 Directogram'})
    elseif i == 2
        title('X03 Directogram')
    else
        title('X04 Directogram')
        xlabel('Date/Time')
    end
    ylabel(cb,'Wave Direction (degrees)')
    ylabel('Frequency (Hz)')
end
        
    


