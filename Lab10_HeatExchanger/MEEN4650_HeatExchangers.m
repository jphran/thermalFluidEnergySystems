% JustinFrancis
% MEEN 4650, TFES
% Dr. Metzger, U of U
% Lab 10, Shell and Tube Heat Exchangers

clc; clear; close all;

%% Data Analysis
clc;
%raw data
rawData = xlsread('HeatExchanger_DataSheet.xlsx');
P_atm = 661; % mmHg
V_dot_cold = rawData(7:end, 1); %gpm
V_dot_hot = rawData(7:end, 2);
T_hi = rawData(7:end, 3); %F, hot side inlet
T_ho = rawData(7:end, 4); % hot side out
T_ci = rawData(7:end, 5); % cold in
T_co = rawData(7:end, 6); %cold out
T_shell = rawData(7:end, 7); %shell
T_amb = rawData(7:end, 8); %amb air

% 1, convert to SI
P_atm = P_atm * 133.322; %Pa
V_dot_cold = V_dot_cold * 0.000063090196; %m3/s
V_dot_hot = V_dot_hot * 0.000063090196; %m3/s
T_hi = (T_hi-32) * 5/9 + 273.15; % K
T_ho = (T_ho-32) * 5/9 + 273.15; % K
T_ci = (T_ci-32) * 5/9 + 273.15; % K
T_co = (T_co-32) * 5/9 + 273.15; % K
T_shell = (T_shell-32) * 5/9 + 273.15; % K
T_amb = (T_amb-32) * 5/9 + 273.15; % K

% 2, calc avg temps
T_cold_avg = 0.5.*(T_ci+T_co);
T_hot_avg = 0.5.*(T_hi+T_ho);

% 3, compute fluid props
[rho_cold, Cp_cold] = WaterProperties(T_cold_avg); %[kg/m3, J/kgK]
[rho_hot , Cp_hot] = WaterProperties(T_hot_avg);

% 4, calc mass flow rate
m_dot_cold = rho_cold.*V_dot_cold; %kg/s
m_dot_hot = rho_hot.*V_dot_hot;

% 5, calc heat capacities
C_cold = m_dot_cold.*Cp_cold;
C_hot = m_dot_hot.*Cp_hot;
C_min = min([C_cold, C_hot]);
C_max = max([C_cold, C_hot]);
C_r = C_min/C_max;

% 6, calc the heat trans rates
q_h = C_hot.*(T_hi-T_ho);
q_c = C_cold.*(T_ci-T_co);

% 7, calc effectiveness of hx *NOTE: q_h yields a better result
effectiveness = q_h./(C_min.*(T_hi-T_ci));

% 8, calc log mean temp diff
dT1 = (T_hi-T_co);
dT2 = (T_ho-T_ci);
dT_lm = (dT1-dT2)./log(dT1./dT2);

%9, calc avg heat trans coeff on the inside tube surface area, using q_h
L = 9*2.54e-2; %m
t = 0.028*2.54e-2; %m
d_outter = 2.12*2.54e-2; %m
d_inner = d_outter-t; %m
A_i = 2*pi*d_inner/2*L + 2*pi*(d_inner/2)^2; %m2
F = 1; %given correction factor
U_i = q_h./(A_i.*F.*dT_lm);

% 10, calc NTU
NTU = (U_i.*A_i)./C_min;

% 11, calc theoretical effectivness
if (C_r < 1)
    effectiveness_theory = (1 - exp(-NTU.*(1-C_r)))./(1 - C_r.*exp(-NTU.*(1-C_r)));
else 
    effectiveness_theory = (NTU)./(1+NTU);
end

%% 12, est rate of heat trans from shell to surrounding
    % finding radiation heat transfer
sig = 5.6703e-8; %boltzman const, W/m2K4
eps = 0.95; %emissivity, given
q_rad = eps*sig.*(T_shell.^4-T_amb.^4)*pi*d_outter*L;
    
    % finding convective heat trans
T_f = 0.5.*(T_shell + T_amb);
[rho_air, mu_air, kappa_air, Cp_air] = AirProperties(T_f, P_atm.*ones(length(T_f)));
g = 9.81; %grav const, m/s2
vu_air = mu_air./rho_air; % kine viscosity
alpha_air = kappa_air./(rho_air.*Cp_air); %thermal diffusivity
Beta = 1/T_f; %vol thermal expansion coeff
Ra_D = (g.*Beta.*(T_shell-T_amb).*d_outter^3)./(vu_air.*alpha_air);
Nu_bar = 0.48.*Ra_D.^(1/4);
h_bar = kappa_air.*Nu_bar/d_outter;
q_conv = h_bar.*pi*d_outter*L.*(T_shell-T_amb);
    
% 13, est percent uncertainty in measured heat trans rates
dT_cold = T_co - T_ci;
dT_hot = T_hi - T_ho;
sig_dT = 0.1;
sig_V_dot = 0.2 * 0.000063090196; % gpm -> m3/s
percent_uncert_cold = ((sig_V_dot./V_dot_cold).^2 + (sig_dT./dT_cold).^2).^(1/2).*100; % percent uncertainty cold stream
percent_uncert_hot = ((sig_V_dot./V_dot_hot).^2 + (sig_dT./dT_hot).^2).^(1/2).*100; % percent uncertainty cold stream

    
    
    
    