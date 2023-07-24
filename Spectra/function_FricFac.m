function[fw,u_br,TED_N] = function_FricFac(See,ff,h,kw,a1,a2,a3,Method)
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
%       - u_br: the representative bottom orbital velocity
%       - TED_N: the energy dissipation(theoretical) using the Neilson
%                 friction factors
%
% Created By: Noah Clark       Last Updated: 7/24/2023
% 
%


%% Preliminaries
loop = size(See);
df = ff(2) - ff(1);
omegaf = 2.*pi.*ff;
rho = 1025;

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
         
                %u_b(j,i) = H(j,i)./(sinh(h(i).*k)*ff(j));  
                u_b(j,i) = (H(j,i)/2)*omegaf(j)/sinh(k*h(i));
            end
        end

        u_bsqr = u_b.^2;
        
        fw = zeros(loop(1),loop(2));
        u_br = zeros(1,loop(2));
        for i = 1:loop(2)
            u_br(i) = sqrt(sum(u_bsqr(:,i)));

            for j = 1:loop(1)
                fw(j,i) = exp(a1.*(u_br(i)./(kw.*omegaf(j))).^a2 + a3);      
            end
        end
        


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
        u_bTp = zeros(1,loop(2));
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
            u_bTp(i) = (Hs(i)/2)*2*pi/(sinh(h(i)*k(i))*Tp);

            for j = 1:loop(1)
                fw(j,i) = exp(a1.*(u_bTp(i)./(kw.*omegaf(j))).^a2 + a3);
            end
        end
        
        u_br = u_bTp;





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
        
        fw = zeros(loop(1),loop(2));
        u_bTm = zeros(1,loop(2));
        for i = 1:loop(2)
                % Solving linear dispersion relationship:
            Lprev = 1; Lnew = 0; thresh = 0.01; delta = 1;
            while delta > thresh
                Lprev = Lnew;
                Lnew = (9.81*Tm(i)^2)/(2*pi)*tanh((2*pi*h(i))/Lprev);
                delta = abs(Lnew - Lprev);
            end
            k = (2*pi)/Lnew;
            %km(i) = function_KwavecalculateSI(Tm(i),h(i));
            u_bTm(i) = ((Hs(i)/2)*2*pi)/(sinh(h(i)*k)*Tm);

            for j = 1:loop(1)
                fw(j,i) = exp(a1.*(u_bTm(i)./(kw.*omegaf(j))).^a2 + a3);
            end
        end
        
        u_br = u_bTm;


end


%% Calculating TED_N
% Calculate the theoretical energy dissipation using the calculated
% friction factors
% Note: Urms is what is origionally used in the final equation but
%        instead we use u_br

% u_b = zeros(loop(1),loop(2));
% for j = 1:loop(1)
%     for i = 1:loop(2)
%         % Solving linear dispersion relationship:
%         Lprev = 1; Lnew = 0; thresh = 0.01; delta = 1;
%         while delta > thresh
%             Lprev = Lnew;
%             Lnew = (9.81*T(j)^2)/(2*pi)*tanh((2*pi*h(i))/Lprev);
%             delta = abs(Lnew - Lprev);
%         end
%         k = (2*pi)/Lnew;
%         
%         u_b(j,i) = (H(j,i)/2)*omegaf(j)/sinh(k*h(i));
%     end
% end

%u_bsqr = u_b.^2;

%u_rms = zeros(1,loop(2));  Urms is what is origionally used in the final
%                           equation but instead we use u_br

TED_N = zeros(loop(1),loop(2));
for i = 1:loop(2)
   % u_rms(i) = sqrt((1/loop(1)).*sum(u_bsqr(:,i)));
    for j = 1:loop(1)
        TED_N(j,i) = 0.6.*rho.*fw(j,i).*u_br(i).^3;
    end
end




end





