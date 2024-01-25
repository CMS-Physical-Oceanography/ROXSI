%% pmBOV.m


    % Noah Clark
    % Created: 6/21/2023
    
    % Purpose: 
    %           - to determine the peak and mean time-weighted frequencies,
    %             periods, bottom orbital velocities, and wavelengths for 
    %             all buoys at Asilomar and China Rock
    %
    
    
%%

clear;clc;

    %Load in data from WBvariables.mat:
load('WBvariables.mat')


%%

           
    % For Asilomar Buoys:
for i = 1:3
    for xx = 1:832
        %Determining L from peak frequency/period:
            [m,indPeakFreq] = max(XSee{i}(:,xx));
            XTp{i}(xx) = 1/Xfreq{i}(indPeakFreq);
            Xfp{i}(xx) = 1/XTp{i}(xx);
            [XpWavelength{i}(xx),XpWaveNumber{i}(xx),XpCelerity{i}(xx)] = ...
            function_wavecalculateSI(XTp{i}(xx),XGivenHsig{i}(xx),Xdepth{i}(xx));
        %Determining L from energy weighted mean frequency/period:
            m0 = trapz(Xfreq{i},XSee{i}(:,xx),1);
            m1 = trapz(Xfreq{i},XSee{i}(:,xx).*Xfreq{1},1);
            XTm{i}(xx) = m0/m1;
            Xfm{i}(xx) = 1/XTm{i}(xx);
            [XmWavelength{i}(xx),XmWaveNumber{i}(xx),XmCelerity{i}(xx)] = ...
            function_wavecalculateSI(XTm{i}(xx),XGivenHsig{i}(xx),Xdepth{i}(xx));
        %Determing peak bottom orbital velocity:
            XpBOV{i}(xx) = (XGivenHsig{i}(xx)*pi)/(XTp{i}(xx)*...
                sinh(XpWaveNumber{i}(xx)*Xdepth{i}(xx)));
        %Determing mean bottom orbital velocity:
            XmBOV{i}(xx) = (XGivenHsig{i}(xx)*pi)/(XTm{i}(xx)*...
                sinh(XmWaveNumber{i}(xx)*Xdepth{i}(xx)));
    end
end



    % For China Rock Buoys:
i = 1;  %Only for buoy B01
for xx = 1:725  %smaller arrays because of time taken out of water
    %Determining L from peak frequency/period:
        [m,indPeakFreq] = max(BSee{i}(:,xx));
        BTp{i}(xx) = 1/Bfreq{i}(indPeakFreq);
        Bfp{i}(xx) = 1/BTp{i}(xx);
        [BpWavelength{i}(xx),BpWaveNumber{i}(xx),BpCelerity{i}(xx)] = ...
            function_wavecalculateSI(BTp{i}(xx),BGivenHsig{i}(xx),Bdepth{i}(xx));
    %Determining L from energy weighted mean frequency/period:
        m0 = trapz(Bfreq{i},BSee{i}(:,xx),1);
        m1 = trapz(Bfreq{i},BSee{i}(:,xx).*Bfreq{1},1);
        BTm{i}(xx) = m0/m1;
        Bfm{i}(xx) = 1/BTm{i}(xx);
        [BmWavelength{i}(xx),BmWaveNumber{i}(xx),BmCelerity{i}(xx)] = ...
            function_wavecalculateSI(BTm{i}(xx),BGivenHsig{i}(xx),Bdepth{i}(xx));
    %Determing peak bottom orbital velocity:
        BpBOV{i}(xx) = (BGivenHsig{i}(xx)*pi)/(BTp{i}(xx)*...
            sinh(BpWaveNumber{i}(xx)*Bdepth{i}(xx)));
    %Determing mean bottom orbital velocity:
        BmBOV{i}(xx) = (BGivenHsig{i}(xx)*pi)/(BTm{i}(xx)*...
            sinh(BmWaveNumber{i}(xx)*Bdepth{i}(xx)));
end

for i = 2:3
    for xx = 1:832
        %Determining L from peak frequency/period:
            [m,indPeakFreq] = max(BSee{i}(:,xx));
            BTp{i}(xx) = 1/Bfreq{i}(indPeakFreq);
            Bfp{i}(xx) = 1/BTp{i}(xx);
            [BpWavelength{i}(xx),BpWaveNumber{i}(xx),BpCelerity{i}(xx)] = ...
            function_wavecalculateSI(BTp{i}(xx),BGivenHsig{i}(xx),Bdepth{i}(xx));
        %Determining L from energy weighted mean frequency/period:
            m0 = trapz(Bfreq{i},BSee{i}(:,xx),1);
            m1 = trapz(Bfreq{i},BSee{i}(:,xx).*Bfreq{1},1);
            BTm{i}(xx) = m0/m1;
            Bfm{i}(xx) = 1/BTm{i}(xx);
            [BmWavelength{i}(xx),BmWaveNumber{i}(xx),BmCelerity{i}(xx)] = ...
            function_wavecalculateSI(BTm{i}(xx),BGivenHsig{i}(xx),Bdepth{i}(xx));
        %Determing peak bottom orbital velocity:
            BpBOV{i}(xx) = (BGivenHsig{i}(xx)*pi)/(BTp{i}(xx)*...
                sinh(BpWaveNumber{i}(xx)*Bdepth{i}(xx)));
        %Determing mean bottom orbital velocity:
            BmBOV{i}(xx) = (BGivenHsig{i}(xx)*pi)/(BTm{i}(xx)*...
                sinh(BmWaveNumber{i}(xx)*Bdepth{i}(xx)));
    end
end
XTp{2}(173) = NaN;  %set to NaN; otherwise would be infinity and mess other things up for the future
Xfp{2}(173) = NaN;

% - China Rock buoy B01 doesn't work for this method because there are only
%    725 values for many of the variables (I'm not sure how to get trapz to 
%    work for this)



%%

    % Saving Variables:
save('WBvariables.mat','-append','BTp','Bfp','BpWavelength',...
    'BpCelerity','BTm','Bfm','BmWavelength','BmCelerity','BpBOV',...
    'BmBOV','XTp','Xfp','XpWavelength','XpCelerity','XTm','Xfm',...
    'XmWavelength','XmCelerity','XpBOV','XmBOV')




