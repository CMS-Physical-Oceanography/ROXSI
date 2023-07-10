function[fw] = function_FricFac(See,ff,h,kw,a1,a2,a3,Method)
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


%%
%function[f_w] = function_FricFac2(See,ff,h,kw,[a1,a2,a3])

% This will eventually become that function
%       - it's not for now just for testing

%% (this should be DELETED once function completed)
% clc;clear;
% load('WBvariables.mat')

%% For now, assume:
% Would be entered into equation:
% See = XSee{1};
% h = Xdepth{1};
% df = 0.0098;
% kw = 0.16;  %roughness length
% ff = (1:129).*df;   %vector of frequencies contained in See
% a1 = 5.5; a2 = -0.2; a3 = -6.3; %Nielson's constants
% - the method should also be an input
%%









%% Preliminaries
loops = size(See);
df = ff(2) - ff(1);
omegaf = 2*pi*ff;

H = 4.*sqrt(See.*df); % CHECK df GOES OUTSIDE SQRT or if there's a *4 (I THINK WHAT I HAVE IS RIGHT)
T = 1./ff;

for i = 1:loops(2)
    Hs(i) = 4.*sqrt(sum(See(:,i)).*df);
end   

switch Method
    
    
%% Representative Orbital Velocity Method

    case {'R','r'}     
        for j = 1:loops(1)
            for i = 1:loops(2)
                k(j,i) = function_KwavecalculateSI(T(j),h(i));  %HAVE A DIFFERENT WAY SO USERS DON'T HAVE TO DOWNLOAD TWO FUNCTIONS
                u_b(j,i) = H(j,i)./(sinh(h(i).*k(j,i)).*ff(j));   
            end
        end

        u_bsqr = u_b.^2;

        for i = 1:loops(2)
            u_br(i) = sqrt(sum(u_bsqr(:,i)));

            for j = 1:loops(1)
                fw(j,i) = exp(a1.*(u_br(i)./(kw.*omegaf(j))).^a2 + a3);      
            end
        end


        % max(u_br);

        %       FOR CALCULATING U_BR:
        % - Originally thought: Times 4 and df inside: max(u_br) = 55
        % - Maybe: Times 4 and df outside: max(u_br) = 5.37
        % - Grimes method: NOT times 4 and df inside : max(u_br) = 13.56

        %       Example Plot (for representative method):
        % - Small energy @ time 263
        % - Medium Energy @ time 189
        % - Large energy @ time 73
        figure(1);clf;
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
        
        for i = 1:loops(2)
            [m,indPeakFreq] = max(See(:,i));
            Tp(i) = 1/ff(indPeakFreq);
            fp(i) = 1/Tp(i);
        end

        for i = 1:loops(2)
            kp(i) = function_KwavecalculateSI(Tp(i),h(i));
            u_bTp(i) = Hs(i)/(sinh(h(i)*kp(i))*fp(i));

            for j = 1:loops(1)
                fw(j,i) = exp(a1.*(u_bTp(i)./(kw.*omegaf(j))).^a2 + a3);
            end
        end

        %   Example Plot (for peak method):
        figure(2);clf;
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
        
 
        for i = 1:loops(2)
            m0 = trapz(ff,See(:,i),1);
            m1 = trapz(ff,See(:,i).*ff',1);
            Tm(i) = m0/m1;
            fm(i) = 1/Tm(i);
        end

        for i = 1:loops(2)
            km(i) = function_KwavecalculateSI(Tm(i),h(i));
            u_bTm(i) = Hs(i)/(sinh(h(i)*km(i))*fm(i));

            for j = 1:loops(1)
                fw(j,i) = exp(a1.*(u_bTm(i)./(kw.*omegaf(j))).^a2 + a3);
            end
        end

        %   Example Plot (for mean method):
        figure(3);clf;
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





