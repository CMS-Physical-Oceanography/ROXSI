function[fw] = function_FricFac(See,ff,h,kw,a1,a2,a3,Method)
%
%FUNCTION_FRICFAC determines the friction factor at each time and at each
% frequency using Neilson's [1992] method taken from "Spectral wave
% dissipation over a barrier reef". A method must be selected which
% specifies which orbital velocity will be calculated.
% 
% Inputs:
%       - See: S eta eta matrix
%       - ff: a vector of frequencies which See follows 
%       - h: water depth vector as a function of time
%       - kw: the roughness length (kw = 0.16 used in paper)
%       - a1, a2, & a3: Neilson's [1992] coefficients which are typically
%         a1 = 5.5, a2 = -0.2, and a3 = -6.3
%       - Method: specify which bottom orbital velocity will be used to
%          perform the calculation
%               - The 3 Options:
%                       - "R"/"r": representative bottom orbital velocity
%                           method (representative used in the paper)
%                       - "M"/"m": mean energy-weighted bottom orbital
%                           velocity method
%                       - "P"/"p": peak bottom orbital velocity method
% 
% Outputs:
%       - fw: the friction factor at each time and at each given frequency
% 
%
% Created By: Noah Clark       Last Updated: 7/11/2023
% 
%
% Still TO-DO:
%               - better comments?
%               - define variables?
%

%% Preliminaries
loop = size(See);
df = ff(2) - ff(1);
omegaf = 2*pi*ff;

H = 4.*sqrt(See.*df); 
T = 1./ff;

Hs = zeros(1,loop(2));
for i = 1:loop(2)
    Hs(i) = 4.*sqrt(sum(See(:,i)).*df);
end   

switch Method
    
    
%% Representative Orbital Velocity Method

    case {'R','r'}     
        u_b = zeros(loop(1),loop(2));
        for j = 1:loop(1)
            for i = 1:loop(2)
                    % Solving linear dispersion relationship:
                Lprev = 1; Lnew = 0; thresh = 0.01; delta = 1;
                while delta > thresh
                    Lprev = Lnew;
                    Lnew = (9.81*T(j)^2)/(2*pi)*tanh((2*pi*h(i))/Lprev);
                    delta = abs(Lnew - Lprev);
                end
                k = (2*pi)/Lnew;
                
                %k(j,i) = function_KwavecalculateSI(T(j),h(i));  %HAVE A DIFFERENT WAY SO USERS DON'T HAVE TO DOWNLOAD TWO FUNCTIONS
                u_b(j,i) = H(j,i)./(sinh(h(i).*k).*ff(j));   
            end
        end

        u_bsqr = u_b.^2;
        
        fw = zeros(loop(1),loop(2));
        for i = 1:loop(2)
            u_br = sqrt(sum(u_bsqr(:,i)));

            for j = 1:loop(1)
                fw(j,i) = exp(a1.*(u_br./(kw.*omegaf(j))).^a2 + a3);      
            end
        end


        %       Example Plot (for representative method):
        % - Small energy @ time 263
        % - Medium Energy @ time 189
        % - Large energy @ time 73
        figure;clf;
        plot(ff,fw(:,263),'g','LineWidth',2)
        hold on
        plot(ff,fw(:,189),'b','LineWidth',2)
        plot(ff,fw(:,73),'r','LineWidth',2)
        xlabel('Frequency (Hz)');ylabel('Friction Factor (fw)');
        legend('Low Energy Time','Medium Energy Time',...
            'High Energy Time','location','northwest')
        title('Representative Orbital Velocity Method')
        grid on



%% Peak Bottom Orbital Velocity Method

    case {'P','p'}
        Tp = zeros(1,loop(2));
        fp = zeros(1,loop(2));
        for i = 1:loop(2)
            [~,indPeakFreq] = max(See(:,i));
            Tp(i) = 1/ff(indPeakFreq);
            fp(i) = 1/Tp(i);
        end
        
        fw = zeros(loop(1),loop(2));
        for i = 1:loop(2)
                % Solving linear dispersion relationship:
            Lprev = 1; Lnew = 0; thresh = 0.01; delta = 1;
            while delta > thresh
                Lprev = Lnew;
                Lnew = (9.81*Tp(i)^2)/(2*pi)*tanh((2*pi*h(i))/Lprev);
                delta = abs(Lnew - Lprev);
            end
            k = (2*pi)/Lnew;
            %kp(i) = function_KwavecalculateSI(Tp(i),h(i));
            u_bTp = Hs(i)/(sinh(h(i)*k(i))*fp(i));

            for j = 1:loop(1)
                fw(j,i) = exp(a1.*(u_bTp./(kw.*omegaf(j))).^a2 + a3);
            end
        end

        %   Example Plot (for peak method):
        figure;clf;
        plot(ff,fw(:,263),'g','LineWidth',2)
        hold on
        plot(ff,fw(:,189),'b','LineWidth',2)
        plot(ff,fw(:,73),'r','LineWidth',2)
        xlabel('Frequency (Hz)');ylabel('Friction Factor (fw)');
        legend('Low Energy Time','Medium Energy Time',...
            'High Energy Time','location','northwest')
        title('Peak Orbital Velocity Method')
        grid on



%% Mean Energy-Weighted Bottom Orbital Velocity Method:  

    case {'M','m'}
        
        Tm = zeros(1,loop(2));
        fm = zeros(1,loop(2));
        for i = 1:loop(2)
            m0 = trapz(ff,See(:,i),1);
            m1 = trapz(ff,See(:,i).*ff',1);
            Tm(i) = m0/m1;
            fm(i) = 1/Tm(i);
        end
        
        km = zeros(1,loop(2));
        u_bTm = zeros(1,loop(2));
        fw = zeros(loop(1),loop(2));
        for i = 1:loop(2)
                % Solving linear dispersion relationship:
            Lprev = 1; Lnew = 0; thresh = 0.01; delta = 1;
            while delta > thresh
                Lprev = Lnew;
                Lnew = (9.81*T(i)^2)/(2*pi)*tanh((2*pi*h(i))/Lprev);
                delta = abs(Lnew - Lprev);
            end
            k = (2*pi)/Lnew;
            %km(i) = function_KwavecalculateSI(Tm(i),h(i));
            u_bTm = Hs(i)/(sinh(h(i)*k)*fm(i));

            for j = 1:loop(1)
                fw(j,i) = exp(a1.*(u_bTm./(kw.*omegaf(j))).^a2 + a3);
            end
        end

        %   Example Plot (for mean method):
        figure;clf;
        plot(ff,fw(:,263),'g','LineWidth',2)
        hold on
        plot(ff,fw(:,189),'b','LineWidth',2)
        plot(ff,fw(:,73),'r','LineWidth',2)
        xlabel('Frequency (Hz)');ylabel('Friction Factor (fw)');
        legend('Low Energy Time','Medium Energy Time',...
            'High Energy Time','location','northwest')
        title('Energy-Weighted Mean Orbital Velocity Method')
        grid on

end


end





