function[f_w] = function_FricFac(Freq,Buoy,Method)
%FUNCTION_FRICFAC determines the friction factor based on a specified
% frequency, buoy, and method; and a plot is created of the friction factor
% curve as a function of frequency. 
%
% Notes:
%       - if 0 is entered as the input frequency, the returned plot will
%          not have a specific frequency and/or friction factor marked
%
%
% Inputs:
%        - Freq: frequency (Hz)
%        - Buoy: the specificed buoy (Ex: 'X01')
%        - Method: the method that should be used to determine the friction
%                  factor (Options: representative, mean, and peak)
% Outputs:
%        - fw: friction factor
%
%

% Loading in variables from WBvariables.mat
load('WBvariables.mat')


% Constants:
A_omegaf = 2*pi*Freq;
kw = 0.16; %m      %horizontal roughness length    %set to 1 in Lowe paper
df = 0.0098;
ff = [1:129].*df;
omegaf = 2*pi*ff;
    % Nielsen(1992) constants:
    a1 = 5.5;
    a2 = -0.2;
    a3 = -6.3;

    
    % Select the buoy to calculate for:
switch Buoy
%     case 'B01'        % at the time, won't work for B01
%         See = BSee{1};
%         h = BavgD{1};
%         Hs = BGivenHsig{1};
    case 'B03'
        See = BSee{2};
        h = BavgD{2};
        Hs = nanmean(BGivenHsig{2});
        Tm = nanmean(BTm{2});
        Tp = nanmean(BTp{2});
    case 'B05'
        See = BSee{3};
        h = BavgD{3};
        Hs = nanmean(BGivenHsig{3});
        Tm = nanmean(BTm{3});
        Tp = nanmean(BTp{3});
    case 'X01'
        See = XSee{1};
        h = XavgD{1};
        Hs = nanmean(XGivenHsig{1});
        Tm = nanmean(XTm{1});
        Tp = nanmean(XTp{1});
    case 'X03'
        See = XSee{2};
        h = XavgD{2};
        Hs = nanmean(XGivenHsig{2});
        Tm = nanmean(XTm{2});
        Tp = nanmean(XTp{2});
    case 'X04'
        See = XSee{3};
        h = XavgD{3};
        Hs = nanmean(XGivenHsig{3});
        Tm = nanmean(XTm{3});
        Tp = nanmean(XTp{3});
end

    
    % Select the method to use to calculate: 
switch Method
    case {'Representative','representative'}
        for j = 1:129
            f(j) = j*df;
            T(j) = 1/f(j);
            omega(j) = 2*pi*f(j);
            H(j) = 4*nansum(See(j,1:832))*df;    %times 4???
            a(j) = H(j)/2;
            [L,k(j),WDP,WS,RD,C] = function_wavecalculateSI(T(j),H(j),h);
            u_bj(j) = a(j)*omega(j)/(sinh(k(j)*h)); %eqn 22 Lowe
        end
        u_bjsq = u_bj.^2;
        u_br = sqrt(sum(u_bjsq));
        f_w = exp(a1.*(u_br./(kw*omegaf)).^a2 + a3);
        
        figure(200);clf;
        plot(ff,f_w)
        title({'Friction Factor vs Frequency using the Representative',...
            sprintf('Bottom Orbital Velocity Method for Buoy %1s',Buoy)},...
            'fontsize',11)
        %xlim([0 0.5]);ylim([0 1]);
        hold on
        
        if A_omegaf ~= 0
            Af_w = exp(a1.*(u_br./(kw*A_omegaf)).^a2 + a3);
        end

    case {'Mean','mean'}
        [L,k,WDP,WS,RD,C] = function_wavecalculateSI(Tm,Hs,h);
        fm = 1/Tm;
        u_bTm = Hs/(sinh(h*k));
        f_w = exp(a1.*(u_bTm./(kw*omegaf)).^a2 + a3);
        
        figure(200);clf;
        plot(ff,f_w)
        title({'Friction Factor vs Frequency using the Mean',...
            sprintf('Bottom Orbital Velocity Method for Buoy %1s',Buoy)},...
            'fontsize',11)
        %xlim([0 0.5]);ylim([0 30]);
        hold on
        
        if A_omegaf ~= 0
            Af_w = exp(a1.*(u_bTm./(kw*A_omegaf)).^a2 + a3);
        end
        
    case {'Peak','peak'}
        [L,k,WDP,WS,RD,C] = function_wavecalculateSI(Tp,Hs,h);
        fp = 1/Tp;
        u_bTp = Hs/(sinh(h*k));
        f_w = exp(a1.*(u_bTp./(kw*omegaf)).^a2 + a3);
        
        figure(200);clf;
        plot(ff,f_w)
        title({'Friction Factor vs Frequency using the Peak',...
            sprintf('Bottom Orbital Velocity Method for Buoy %1s',Buoy)},...
            'fontsize',11)
        %xlim([0 0.5]);ylim([0 3]);
        hold on
        
        if A_omegaf ~= 0
            Af_w = exp(a1.*(u_bTp./(kw*A_omegaf)).^a2 + a3);
        end
end

grid on
xlabel('Frequency, {\it f} (Hz)')
ylabel('Friction Factor {\it f_w}')

if A_omegaf ~= 0
    plot(Freq,Af_w,'.m','Markersize',15)
    legend('Friction Factor Curve','Frequency of Interest','location','best')

    ann1 = annotation('textbox',[0.14,0.8,0.37,0.07],'string',...
        sprintf('Selected Frequency, f (Hz): %g\n',round(Freq,3)),'fontsize',10);
    ann2 = annotation('textbox',[0.14,0.73,0.37,0.07],'string',...
        sprintf('Friction Factor, f_w: %g\n',round(Af_w,3)),'fontsize',10);
    ann1.BackgroundColor = 'w';
    ann2.BackgroundColor = 'w';

    f_w = Af_w
else
    f_w = 'Not Available'
end

end

% I think something about the wave heights is wrong for mean and peak


