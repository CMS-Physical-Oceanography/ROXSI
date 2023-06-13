%% Spectrograms for Every Buoy
% Noah Clark
% 6/9/2023

clc;clear;

% Load in data from WBdata function:
    [Btime,Bfreq,BSee,BGivenHsig,Bdepth,BEMEM,Xtime,Xfreq,XSee,...
        XGivenHsig,Xdepth,XEMEM] = load_WBdata(1);

% Creating the lines to go on spectrogram to say "No data for this time":
    LB1x1 = [datetime('02-Jul-2022 00:00:00','TimeZone','America/Los_Angeles')...
        , datetime('02-Jul-2022 00:00:00','TimeZone','America/Los_Angeles')];
    LB1x2 = [datetime('06-Jul-2022 12:00:00','TimeZone','America/Los_Angeles')...
        , datetime('06-Jul-2022 12:00:00','TimeZone','America/Los_Angeles')];
    LB1y1 = [0 0.55];
    LB1y2 = [0 0.55];


% Creating a vector of frequencies from 0.01 Hz to 1.29 Hz:
    Gfreq = (1:129).*0.01;

% Making the spectrograms for the China Rock Buoys:
    for i = 1:3
        figure(i);clf;
        delete(findall(gcf,'type','annotation'))
        set(gcf,'position',[0,570,450,450])
        pcolor(Btime{i},Gfreq,BSee{i})
        shading flat
        cb = colorbar;
        caxis([0 3])
        ylim([0 0.4])
        if i == 1
            hold on
            plot(LB1x1,LB1y1,'-r','LineWidth',2)
            plot(LB1x2,LB1y2,'-r','LineWidth',2)
            title('B01 Spectrogram')
            an = annotation('textbox',[0.4177,0.7,0.08,0.187],'string',...
                'No Data for this Time','fontsize',8);
            set(an,'color','r')
        elseif i == 2
            title('B03 Spectrogram')
        else
            title('B05 Spectrogram')
        end
        xlabel('Date/Time')
        ylabel(cb,'Wave Energy (m^2/Hz)')
        ylabel('Frequency (Hz)')
    end
    
% Making the spectrograms for the Asilomar Buoys:
    for i = 1:3
        figure(i + 3);clf;
        set(gcf,'position',[0,570,450,450])
        pcolor(Xtime{i},Gfreq,XSee{i})
        shading flat
        cb = colorbar;
        caxis([0 3])
        ylim([0 0.4])
        if i == 1
            title('X01 Spectrogram')
        elseif i == 2
            title('X03 Spectrogram')
        else
            title('X04 Spectrogram')
        end
        xlabel('Date/Time')
        ylabel(cb,'Wave Energy (m^2/Hz)')
        ylabel('Frequency (Hz)')
    end
        