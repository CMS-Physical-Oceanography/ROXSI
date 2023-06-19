%% plot_X01animation.m

    % Noah Clark
    % Created: 6/16/2023
    
    % Purpose: 
    %          - create and save the animation of the X01 spectrum 
    %          - determine the significant wave height predicted by the
    %            model using the integration method 
    %          - determining peak period, frequency, wavelength, celerity,
    %            and bottom orbital velocity
    %          - determining energy-weighted mean period, frequency, 
    %            wavelength, celerity, and bottom orbital velocity
    %            
 
    
    
%%

clear;clc;

    %Load in data from WBvariables.mat:
load('WBvariables.mat')


%%
    % Calculate the significant wave height predicted by the model using
    % integration method:
for i = 1:787
    M_Hsig(i) = 4*(sqrt(trapz(ModelSee.X01.f,ModelSee.X01.Snn_f(:,i),1)));
end



        % Creating animation:
    
    %Change the folder to save the animation figures to
cd 'C:\Users\nsc4750\Documents\CMS Summer\Start (5-16)\Start Figures\Youtube Animation'
addpath 'C:\Users\nsc4750\Documents\CMS Summer\Start (5-16)'

figure(1);clf
set(gcf,'position',[1350,40,450,450])
%277:444(June 27th - July 3rd)(Used for Youtube video)
for xx = 1:832
  figure(1);clf;
   delete(findall(gcf,'type','annotation'))
    
    if xx < 215
       annotation('textbox',[0.32,0.8,0.4,0.07],'string','No Wind Data for this Time')
    end
    if xx > 214
       subplot(3,1,1)
       title('Wind Direction and Speed (m/s)')
       hidden_arrow = compass(9.5,0);
       hidden_arrow.Color = 'none';
       hold on
       [u,v] = pol2cart(deg2rad(WindDir(xx - 214)),WindSpd(xx - 214));
       c = compass(u,v,'r');
       c.LineWidth = 3;
       ax = ancestor(c(1),'axes');
       th = findall(ax,'Type','Text');
       set(th,'FontSize',8)
       %c.fontsize = 7;
       view(90,-90)
        annotation('textbox',[0.7,0.83,0.27,0.1],'String',...
            sprintf('Wind Speed (m/s): %1s',...
            sprintf('%g\n',round(WindSpd(xx - 214),2))),'fontsize',9)
        annotation('textbox',[0.7,0.73,0.27,0.1],'String',...
            sprintf('Wind Direction (degrees): %1s',...
            sprintf('%g\n',round(WindDir(xx - 214),2))),'fontsize',9)
        title({'Wind Speed and Direction',''},'fontsize',11)
        hold off
    end

%Plotting:
    subplot(3,1,[2:3])
    plot(Xfreq{1},XSee{1}(:,xx),'k','LineWidth',1.5)

    if xx > 35 && xx < 823
        hold on
        yy = xx - 35;
        plot(ModelSee.X01.f,ModelSee.X01.Snn_f(:,yy),'r','LineWidth',1.5)

        annotation('textbox',[0.55,0.49,0.34,0.06],'String',...
            sprintf('Theoretical H_1_/_3 (m): %1s',...
            sprintf('%g\n',round(M_Hsig(yy),2))),'fontsize',9)
        legend('Observed','Modeled','location','southeast')
    end

    hold on
    xlabel('Frequency (Hz)')
    ylabel('Energy (m^2/Hz)')
    grid on
    axis([0 0.5 0 3.5])
    title({'',sprintf('Frequency Spectrum for Bouy X01 @ %1s',Xtime{1}(xx))})
    if xx < 35 && xx > 823 
        legend('Observed','Modeled','location','southeast')
    end    
    
    %Calculations:
        %Determining L from peak frequency/period:
            [m,indPeakFreq] = max(XSee{1}(:,xx));
            Tp(xx) = 1/Xfreq{1}(indPeakFreq);
            fp(xx) = 1/Tp(xx);
            [pWavelength(xx),pWaveNumber(xx),pCelerity(xx)] = ...
            wavecalculateSI(Tp(xx),XGivenHsig{1}(xx),Xdepth{1}(xx));
        %Determining L from energy weighted mean frequency/period:
            m0 = trapz(Xfreq{1},XSee{1}(:,xx),1);
            m1 = trapz(Xfreq{1},XSee{1}(:,xx).*Xfreq{1},1);
            Tm(xx) = m0/m1;
            fm(xx) = 1/Tm(xx);
            [mWavelength(xx),mWaveNumber(xx),mCelerity(xx)] = ...
            wavecalculateSI(Tm(xx),XGivenHsig{1}(xx),Xdepth{1}(xx));
        %Determing peak bottom orbital velocity:
            pBOV(xx) = (XGivenHsig{1}(xx)*pi)/(Tp(xx)*...
                sinh(pWaveNumber(xx)*Xdepth{1}(xx)));
        %Determing mean bottom orbital velocity:
            mBOV(xx) = (XGivenHsig{1}(xx)*pi)/(Tm(xx)*...
                sinh(mWaveNumber(xx)*Xdepth{1}(xx)));
                
        %Adding Annotations:
            annotation('textbox',[0.55,0.55,0.34,0.06],'String',...
               sprintf(' Measured H_1_/_3 (m): %1s',...
               sprintf('%g\n',round(Xt_Hsig{1}(xx),2))),'fontsize',9)
           annotation('textbox',[0.65,0.41,0.229,0.06],'String',...
               sprintf('L_p (m): %1s',...
               sprintf('%g\n',round(pWavelength(xx),2))),'fontsize',9)
            annotation('textbox',[0.65,0.35,0.229,0.06],'String',...
               sprintf('L_m (m): %1s',...
               sprintf('%g\n',round(mWavelength(xx),2))),'fontsize',9)
           annotation('textbox',[0.65,0.29,0.229,0.06],'String',...
               sprintf('u_p (m/s): %1s',...
               sprintf('%g\n',round(pBOV(xx),2))),'fontsize',9)
           annotation('textbox',[0.65,0.23,0.229,0.06],'String',...
               sprintf('u_m (m/s): %1s',...
              sprintf('%g\n',round(mBOV(xx),2))),'fontsize',9)
    
%To display the progress of the for-loop:
    if xx == 0
        disp('X01 Animation Loop is Starting')
    elseif xx == 83
        disp('X01 Animation Loop is 1/10 Complete')
    elseif xx == 166
        disp('X01 Animation Loop is 2/10 Complete')
    elseif xx == 249
        disp('X01 Animation Loop is 3/10 Complete')
    elseif xx == 332
        disp('X01 Animation Loop is 4/10 Complete')
    elseif xx == 415
        disp('X01 Animation Loop is 5/10 Complete')
    elseif xx == 500
        disp('X01 Animation Loop is 6/10 Complete')
    elseif xx == 580
        disp('X01 Animation Loop is 7/10 Complete')
    elseif xx == 664
        disp('X01 Animation Loop is 8/10 Complete')
    elseif xx == 747
        disp('X01 Animation Loop is 9/10 Complete')
    elseif xx == 832
        disp('X01 Animation Loop is Complete')
    end
    
%%pause(0.3);clf;    %this is for running the animation in matlab
%%saveas(figure(1),sprintf('YT1_Frame%1d.jpeg',xx))
end





