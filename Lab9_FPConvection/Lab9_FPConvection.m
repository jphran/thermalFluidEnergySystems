% Justin Francis
% MEEN 4650, TFES
% Dr. Metzger, U of U
% Lab 9, Flat Plate Convection

clc; clear; close all;

data = load('ConvectionData.dat');
T_s = data(:, 2) + 273.15; %[K]
thermocoupleNum = data(:,1);

%% Calcs from measured data (1)
%plate dims
zeta = 77 * 1e-3; %distance to leading edge[m]
L_h = 153 * 1e-3; %length of plate[m]
L = 230 * 1e-3; %leading edge + plate[m]
w = 68 * 1e-3; %plate[m]
w_T = 144 * 1e-3; %teflon[m]
L_t = (230 + 67) * 1e-3; %teflon[m]
t = 13.9 * 1e-3; %plate thickness[m]
%strip heater dims (mounted side-by-side or width-to-width)
L_sh = 76.2 * 1e-3; %[m]
w_sh = 25.4 * 1e-3; %[m]
R = 156.7; %total resistance of heaters in parrallel[ohm]
V_AC = 45.5; %AC voltage[VAC]
loc = [85 92 102 112 123 123 134 143 153 162 173 173 186 196 209 219] * 1e-3; %[m] from leading edge
%measured dims
fanFreq = 9; %[Hz]
T_amb = 273.15 + 20.8; %[K]
T_inf = T_amb;
P_amb = 87.939447; %[kPa]
P_dyn = 0.0584774; %[kPa]
[rho,mu,k,Cp] = AirProperties(T_amb, P_amb*1e3);
V_fs = sqrt(2*P_dyn/rho); %free stream velocity[m/s]

%% 1, calcs from measured data
%a, heat flux from top surf
netHeatFlux_top = V_AC^2/(2*R*L_h*w); %perfect efficiency[W/m^2]
%b, heat trans coeff
localHeatTransCoeff = netHeatFlux_top./(T_s-T_inf); %[W/(m^2*K)]
avgHeatTransCoeff = 1/(loc(length(loc))-loc(1)) * trapz(localHeatTransCoeff, loc);
%c, nusselt num
T_f = (T_s + T_inf)/2; %film temp[K]
[rho_f, mu_f, k_f, Cp_f] = AirProperties(T_f, P_amb*1e3);
localNusselt = localHeatTransCoeff.*loc./k_f; 
avgT_f = mean(T_f);
[rho_bar, mu_bar, k_bar, Cp_bar] = AirProperties(avgT_f, P_amb*1e3);
avgNusselt = avgHeatTransCoeff*L./k_bar;

%% 2, theoretical predictions
%a, local nusselt and heat 
Re_x = (V_fs.*loc)./(mu_bar/rho_bar);
Pr = (mu_bar/rho_bar)/(k_bar/(Cp_bar*rho_bar));
Re_cr = 5e5;

% theoreticalLocalNusselt_laminar = (0.453*Re_x.^(0.5)*Pr^(1/3))/(1 - (zeta./loc).^(3/4)).^(1/3);
% theoreticalLocalNusselt_turbulent = (0.031*Re_x.^(4/5)*Pr^(1/3))/(1 - (zeta./loc).^(9/10)).^(1/9);
theoreticalLocalNusselt = zeros(length(Re_x),1);

for i = 1:length(Re_x)
    if Re_x(i) <= Re_cr
        theoreticalLocalNusselt(i) = (0.453*Re_x(i)^(0.5)*Pr^(1/3))/(1 - (zeta/loc(i))^(3/4))^(1/3);
    else
        theoreticalLocalNusselt(i) = (0.031*Re_x(i)^(4/5)*Pr^(1/3))/(1 - (zeta/loc(i))^(9/10))^(1/9);
    end
end









