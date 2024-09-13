function [WorkSpace] = Turbine1D(DesignPar,Constants)

% This function calculates through the entire turbine stage based on two input structures 'DesignPar' and 'Constants'.
% The massflow in station 2 & 3 are matched through an 'iteration' on alpha2 and beta3

% For a more detailed explanation see ppt and Excel file equivalent


%% General Calculations

U = DesignPar.RPM/60*2*pi*DesignPar.Rmean;      % Peripheral speed [m/s]



%% STATION 1

A1 = pi*((DesignPar.Rmean+DesignPar.H1/2)^2-(DesignPar.Rmean-DesignPar.H1/2)^2);
P1 = DesignPar.P01/(1+(Constants.Gamma-1)/2*DesignPar.M1^2)^(Constants.Gamma/(Constants.Gamma-1));
T1 = DesignPar.T01/(1+(Constants.Gamma-1)/2*DesignPar.M1^2);
V1 = sqrt(Constants.Gamma*Constants.R*T1)*DesignPar.M1;
Rho1 = P1/Constants.R/T1;
Mu1 = Constants.Sutherlandmu0*(T1/Constants.SutherlandT0)^(3/2)*(Constants.SutherlandT0 + Constants.SutherlandS)/(T1 + Constants.SutherlandS);

Mf1 = Rho1*V1*A1;



%% STATION 2

% 'iteration' needed to match massflow by adjusting alpha2

% Initialize a vector of possible Stator outlet angles
alpha2 = [5:Constants.AngleAccuracy:89];        % Stator exit absolute flow angle (compared to axial) [deg.]

% Calculate all quantities in Station 2
A2 = pi*((DesignPar.Rmean+DesignPar.H1*DesignPar.H2H1/2)^2-(DesignPar.Rmean-DesignPar.H1*DesignPar.H2H1/2)^2); 
P2 = DesignPar.Rp*(DesignPar.P01-DesignPar.P01/DesignPar.PR)+DesignPar.P01/DesignPar.PR;
T02 = DesignPar.T01;
V2is = sqrt(2*Constants.Cp*DesignPar.T01*(1-(P2/DesignPar.P01)^((Constants.Gamma-1)/(Constants.Gamma))));
V2 = sqrt(Constants.Eta_st*V2is^2);
T2 = T02-V2^2/2/Constants.Cp;
M2 = V2/sqrt(Constants.Gamma*Constants.R*T2);
Rho2 = P2/Constants.R/T2;
P02 = P2*(1+(Constants.Gamma-1)/2*M2^2)^(Constants.Gamma/(Constants.Gamma-1));
Mu2 = Constants.Sutherlandmu0*(T2/Constants.SutherlandT0)^(3/2)*(Constants.SutherlandT0 + Constants.SutherlandS)/(T2 + Constants.SutherlandS);
Omega_st = (DesignPar.P01-P02)./(P02-P2);

V2ax = V2.*cos(alpha2./180.*pi);
W2ax = V2ax;
V2t = V2.*sin(alpha2./180.*pi);
W2t =V2t-U;
beta2 = atan2(W2t,W2ax).*180./pi;
W2 = sqrt(W2t.^2+W2ax.^2);
T02r = T2 + W2.^2./2./Constants.Cp;
P02r = P2.*(T02r./T2).^(Constants.Gamma/(Constants.Gamma-1));
M2r = W2./sqrt(Constants.Gamma.*Constants.R.*T2);

% Massflow evaluation
Mf2 = Rho2.*V2ax.*A2;
Error2 = (Mf2-Mf1)./Mf1.*100;                         % error with station 1 massflow [%]

% !! 'Delete' possible second solutions which give rise to supersonic AXIAL Mach numbers !!
for i=1:length(alpha2)
    M2ax = V2ax(i)/sqrt(Constants.Gamma*Constants.R*T2);
    if M2ax > 1
        Error2(i) = 10^10;
    end
end

% Select alpha2 with matching massflow
[minerror minpos] = min(abs(Error2));
alpha2 = alpha2(minpos);

% Attribute the final values
V2ax = V2ax(minpos);
W2ax = W2ax(minpos);
V2t = V2t(minpos);
W2t = W2t(minpos);
beta2 = beta2(minpos);
W2 = W2(minpos);
T02r = T02r(minpos);
P02r = P02r(minpos);
M2r = M2r(minpos);

Mf2 = Mf2(minpos);
Error2 = Error2(minpos);



%% STATION 3

% 'iteration' needed to match massflow by adjusting beta3
beta3 = [-89:Constants.AngleAccuracy:-30];            % Rotor exit relative flow angle (compared to axial) [deg.]

% Calculate all quantities at Station 3
A3 = pi*((DesignPar.Rmean+DesignPar.H1*DesignPar.H2H1*DesignPar.H3H2/2)^2-(DesignPar.Rmean-DesignPar.H1*DesignPar.H2H1*DesignPar.H3H2/2)^2);
T03r = T02r;
P3 = DesignPar.P01/DesignPar.PR;
T3is = T03r./(P02r./P3).^((Constants.Gamma-1)/Constants.Gamma);
W3is = sqrt((T03r-T3is)*Constants.Cp*2);
W3 = sqrt(W3is^2*Constants.Eta_rot);
T3 = T03r-W3.^2./2./Constants.Cp;
M3r = W3./sqrt(Constants.Gamma.*Constants.R.*T3);
P03r = P3.*(1+(Constants.Gamma-1)./2.*M3r.^2).^(Constants.Gamma./(Constants.Gamma-1));
Mu3 = Constants.Sutherlandmu0*(T3/Constants.SutherlandT0)^(3/2)*(Constants.SutherlandT0 + Constants.SutherlandS)/(T3 + Constants.SutherlandS);
Rho3 = P3./Constants.R./T3;
Omega_rot = (P02r-P03r)./(P03r-P3).*100;

W3t = W3.*sin(beta3./180.*pi);
V3t = W3t + U;
V3ax = W3.*cos(beta3./180.*pi);
V3 = sqrt(V3t.^2 + V3ax.^2);
T03 = T3 + V3.^2./Constants.Cp./2;
M3 = V3./sqrt(Constants.Gamma.*Constants.R.*T3);
P03 = P3.*(1+(Constants.Gamma-1)./2.*M3.^2).^(Constants.Gamma./(Constants.Gamma-1));
alpha3 = atan2(V3t,V3ax).*180./pi;

% Massflow evaluation
Mf3 = Rho3.*A3.*V3ax;
Error3 = (Mf1-Mf3)./(Mf1).*100;                          % error with station 1 massflow [%]

% Select angle with matching massflow
[minerror minpos] = min(abs(Error3));

% Attribute the final values
beta3 = beta3(minpos);

W3t = W3t(minpos);
V3t = V3t(minpos);
V3ax = V3ax(minpos);
V3 = V3(minpos);
T03 = T03(minpos);
M3 = M3(minpos);
P03 = P03(minpos);
alpha3 = alpha3(minpos);

Mf3 = Mf3(minpos);
Error3 = Error3(minpos);                    % error with station 1 massflow [%]



%% Overall Quantities

dH = Constants.Cp*(DesignPar.T01-T03);
Power = Mf1*dH;

FlowCoefficient = (V2ax+V3ax)/2/U;
StageLoading = dH/U^2;
RpTrue = (P2-P3)/(P1-P3);
NonDimRotSpeed = U/sqrt(DesignPar.T01);
MeanVelRatio = sqrt(U^2/2/dH);
Efficiency = Power/(Mf1*Constants.Cp*DesignPar.T01*(1-(P03/DesignPar.P01)^((Constants.Gamma-1)/Constants.Gamma)));

ReynoldsStator = Rho2*V2/Mu2;
ReynoldsRotor = Rho3*W3/Mu3;
turning = abs(beta2) + abs(beta3);
if (beta2<0 || beta3>0 || alpha2>79 || beta3>75 || turning>140 || Mf3 >30)
    Efficiency = 0.00001;
    StageLoading = 0.00001;
end

%% Export the full workspace

w = whos;

for a = 1:length(w)
    str.(w(a).name) = eval(w(a).name);
end

WorkSpace = str;

end
