%% plot_WBcalcHsig.m

    % Noah Clark
    % Created: 6/16/2023
    
    % Purpose: 
    %           - use the method where you calculate the area under the
    %             time-averaged spectrum to determine the significant wave
    %             height
    %           - create a table with the 6 average significant wave heights
    %           - create a plot comparing the calculated Hsigs versus the
    %           given Hsigs for the China Rock buoys
 
    
    
%%

clear;clc;

    %Load in data from WBvariables.mat:
load('WBvariables.mat')

df = 0.0098;
BBuoyName = {'B01' 'B03' 'B05'};
XBuoyName = {'X01' 'X03' 'X04'};


%%

    %China Rock (B) - for each hour:
for i = 1:3
    Bt_int{i} = sum(BSee{i}).*df; %the integration
    Bt_Hsig{i} = 4*sqrt(Bt_int{i}); %m
end
    %China Rock (B) - for each buoy:
for i = 1:3
    Bbuoy_int{i} = sum(Bmeanspec{i}).*df; %the integration
    Bbuoy_Hsig{i} = 4*sqrt(Bbuoy_int{i}); %m
end

    %Asilomar (X) - for each hour:
for i = 1:3
    Xt_int{i} = sum(XSee{i}).*df; %the integration
    Xt_Hsig{i} = 4*sqrt(Xt_int{i}); %m
end
    %Asilomar (X) - for each buoy:
for i = 1:3
    Xbuoy_int{i} = sum(Xmeanspec{i}).*df; %the integration
    Xbuoy_Hsig{i} = 4*sqrt(Xbuoy_int{i}); %m
end

    %Create table with the 6 significant wave heights:
Hsig = [Bbuoy_Hsig{1};Bbuoy_Hsig{2};Bbuoy_Hsig{3};Xbuoy_Hsig{1};...
    Xbuoy_Hsig{2};Xbuoy_Hsig{3}];
Buoy = ['B01';'B03';'B05';'X01';'X03';'X04'];
Table_Hsig = table;
Table_Hsig.Buoy = Buoy;
Table_Hsig.Hsig = Hsig

    %Create figure comparing the calculated Hsigs versus the given Hsigs
        %China Rock (B):
        
figure(1);clf;
set(gcf,'position',[900,40,450,450])
x1 = datetime('15-Jun-2022','TimeZone','America/Los_Angeles');
x2 = datetime('21-Jul-2022','TimeZone','America/Los_Angeles');

for i = 1:3
    t = Btime{i};
    subplot(3,1,i)
    %xlim([x1 x2])
    
    plot(t(1:end-1),BGivenHsig{i}(1:end-1),'-b','LineWidth',1.5)
    hold on
    plot(t,Bt_Hsig{i},'-r','LineWidth',1)
    hold on
    if i == 1
        plot(t(397:398),BGivenHsig{i}(397:398),'-w','LineWidth',1.5)
        hold on
        plot(t(397:398),Bt_Hsig{i}(397:398),'-w','LineWidth',1)
        hold on
    end
    ylim([0 4]) 
    xlim([x1,x2])
    xlabel('Time (hrs)')
    ylabel('Hsig (m)')
    title({'',...
        sprintf('Given Versus Calculated Significant Wave Height for %s'...
        ,BBuoyName{i})})
    legend('Given Hsig','Calculated Hsig','location','north','fontsize',6)
    grid on
end


%%

    % Saving Variables:
save('WBvariables.mat','-append','Bt_Hsig','Xt_Hsig')


