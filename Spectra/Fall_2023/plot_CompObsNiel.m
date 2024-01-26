%% plot_CompObsNiel.m
%
%   Noah Clark      11/29/23
%
% Purpose: To plot the integrated obs diss vs the integrated nielsen diss.
%           Integrate over the sea and swell bands separately. 

clc;clear;

%% Load Data

load('TED.mat')

%% Integrate and Assign 

% <0.8 Hz is the Swell Band (ind: <= 8)
% >0.8 Hz is the Sea Band (ind: >8)

for i = 1:8
    Swell_intObsTED{i} = trapz(ff(1:8),Avg_obsTED{i}(1:8));
    Sea_intObsTED{i} = trapz(ff(9:end),Avg_obsTED{i}(9:end));

    Swell_intNielsenTED{i} = trapz(ff(1:8),Avg_NielsonTED{i}(1:8));
    Sea_intNielsenTED{i} = trapz(ff(9:end),Avg_NielsonTED{i}(9:end));

    % Swell_intObsTED{i} = trapz(Avg_obsTED{i}(1:8));
    % Sea_intObsTED{i} = trapz(Avg_obsTED{i}(9:end));
    % 
    % Swell_intNielsenTED{i} = trapz(Avg_NielsonTED{i}(1:8));
    % Sea_intNielsenTED{i} = trapz(Avg_NielsonTED{i}(9:end));
end


% Swell_intObsTED{1} = trapz(ff(1:8),Avg_obsTED{5}(1:8));
% Swell_intObsTED{2} = trapz(ff(1:8),Avg_obsTED{6}(1:8));
% Swell_intObsTED{3} = trapz(ff(1:8),Avg_obsTED{7}(1:8));
% Swell_intObsTED{4} = trapz(ff(1:8),Avg_obsTED{8}(1:8));
% 
% Sea_intObsTED{1} = trapz(ff(9:end),Avg_obsTED{5}(9:end));
% Sea_intObsTED{2} = trapz(ff(9:end),Avg_obsTED{6}(9:end));
% Sea_intObsTED{3} = trapz(ff(9:end),Avg_obsTED{7}(9:end));
% Sea_intObsTED{4} = trapz(ff(9:end),Avg_obsTED{8}(9:end));
% 
% 
% Swell_intNielsenTED{1} = trapz(ff(1:8),Avg_NielsonTED{5}(1:8));
% Swell_intNielsenTED{2} = trapz(ff(1:8),Avg_NielsonTED{6}(1:8));
% Swell_intNielsenTED{3} = trapz(ff(1:8),Avg_NielsonTED{7}(1:8));
% Swell_intNielsenTED{4} = trapz(ff(1:8),Avg_NielsonTED{8}(1:8));
% 
% Sea_intNielsenTED{1} = trapz(ff(9:end),Avg_NielsonTED{5}(9:end));
% Sea_intNielsenTED{2} = trapz(ff(9:end),Avg_NielsonTED{6}(9:end));
% Sea_intNielsenTED{3} = trapz(ff(9:end),Avg_NielsonTED{7}(9:end));
% Sea_intNielsenTED{4} = trapz(ff(9:end),Avg_NielsonTED{8}(9:end));


%% Plot

figure(1);clf;
for i = 1:8
    plot(Sea_intObsTED{i},Sea_intNielsenTED{i},'xr','MarkerSize',10,'LineWidth',3)
    hold on
    plot(Swell_intObsTED{i},Swell_intNielsenTED{i},'xb','MarkerSize',10,'LineWidth',3)
end
grid on
xlabel('Integrated Obs TED')
ylabel('Integrated Nielsen TED')

legend('Sea Band','Swell Band','location','SouthEast')
plot([-0.05 0.5],[-0.05 0.5],'--k','linewidth',1.5)
% xlim([-0.1 1])
% ylim([-0.05 ])
% xlim([-3 48])
% ylim([-5 20])

