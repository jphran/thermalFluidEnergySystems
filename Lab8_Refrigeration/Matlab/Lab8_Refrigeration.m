%Justin Francis
%TFES, lab 8 Refrig
%march 23, 2020

clc; clear; close all;

% load the CoolProp package into Matlab
import py.CoolProp.CoolProp.PropsSI
% set the reference state to match that in the available P-h diagram
py.CoolProp.CoolProp.set_reference_state('R134a','IIR')

Tamb = 23.5 +273.15; %K
Pamb = 651.4 * 0.133322; %kPa
Wfan = 130; %W
flowRate = [.2, .15, .1, .07]; %gpm
totalPower = [840, 680, 570, 515]; %W
T1 = [62, 63, 69, 70]; %F
T2 = [120, 140, 135, 132]; %F
T3 = [98, 82, 72, 75]; %F
T4 = [39, 24, 2, -5]; %F
T5 = [62, 68, 70, 70]; %F
P1 = [31, 20, 10, 5]; %psig
P2 = [147, 135, 136, 108]; %psig
P3 = [145, 135, 114, 105]; %psig
P4 = [40, 27, 13, 8]; %psig
P5 = [35, 23, 12, 9]; %psig
Tc1 = [35.3, 32.8, 29.5, 28.2]; %C
Tc2 = [34.3, 29.7, 25.7, 25]; %C
Te1 = [15.7, 11.7, 11.7, 12.2]; %C
Te2 = [18.9, 22.5, 23, 23]; %C

for i = 1:length(flowRate)
    T1(i) = (T1(i)+459.67) * 5/9; %K
    T2(i) = (T2(i)+459.67) * 5/9; %K
    T3(i) = (T3(i)+459.67) * 5/9; %K
    T4(i) = (T4(i)+459.67) * 5/9; %K
    T5(i) = (T5(i)+459.67) * 5/9; %K
    Tc1(i) = Tc1(i) + 273.15; %K
    Tc2(i) = Tc2(i) + 273.15; %K
    Te1(i) = Te1(i) + 273.15; %K
    Te2(i) = Te2(i) + 273.15; %K
    P1(i) = (P1(i) * 6.89476) + Pamb; %kPa
    P2(i) = (P2(i) * 6.89476) + Pamb; %kPa
    P3(i) = (P3(i) * 6.89476) + Pamb; %kPa
    P4(i) = (P4(i) * 6.89476) + Pamb; %kPa
    P5(i) = (P5(i) * 6.89476) + Pamb; %kPa
    flowRate(i) = flowRate(i) * 0.000063090196; %m3/s
    %crunch
    Te_avg(i) = 0.5 * (Te1(i) + Te2(i)); %K
    Tc_avg(i) = 0.5 * (Tc1(i) + Tc2(i)); %K
end



%% data anal
for i = 1:4
    %3
    H1(i)= py.CoolProp.CoolProp.PropsSI('H','T',T1(i),'P',P1(i)*10^3,'R134a'); %J/kg
    H2(i)= py.CoolProp.CoolProp.PropsSI('H','T',T2(i),'P',P2(i)*10^3,'R134a'); %J/kg
    H3(i)= py.CoolProp.CoolProp.PropsSI('H','T',T3(i),'P',P3(i)*10^3,'R134a'); %J/kg
    H4(i)= H3(i); %py.CoolProp.CoolProp.PropsSI('H','T',T4(i),'P',P4(i),'R134a'); %J/kg
    H5(i)= py.CoolProp.CoolProp.PropsSI('H','T',T5(i),'P',P5(i)*10^3,'R134a'); %J/kg
    %4
    r1(i)=py.CoolProp.CoolProp.PropsSI('D','T',T1(i),'P',P1(i)*10^3,'R134a'); %kg/m3
    r2(i)=py.CoolProp.CoolProp.PropsSI('D','T',T2(i),'P',P2(i)*10^3,'R134a');
    r3(i)=py.CoolProp.CoolProp.PropsSI('D','T',T3(i),'P',P3(i)*10^3,'R134a');
    r4(i)=py.CoolProp.CoolProp.PropsSI('D','T',T3(i),'P',P4(i)*10^3,'R134a');
    r5(i)=py.CoolProp.CoolProp.PropsSI('D','T',T4(i),'P',P5(i)*10^3,'R134a');

end

%5
mdot = r3.*flowRate; %kg/s
%6
Wcomp = totalPower - Wfan; %W
%7
eta = .78; %given
Win = eta.*Wcomp; % W
%8
win = (Win./mdot) .* 10^-3; %kJ/kg
%9
qH = -(H3-H2).*10^-3; %kJ/kg
qL = (H1-H4).*10^-3;
%10
qloss = win + ((H1-H2).*10^-3); %kJ/kg
%11
COPR = qL./win;
%12
for i = 1:4
    s1(i) = py.CoolProp.CoolProp.PropsSI('S', 'T', T1(i), 'P', P1(i)*10^3, 'R134a');
    h2s(i) = py.CoolProp.CoolProp.PropsSI('H', 'S', s1(i), 'T', T2(i), 'R134a');
end
%13
etaC = ((h2s-H1).*10^-3)./(((H2-H1).*10^-3)+qloss) * 100;


%% figures 1a
figure(1);
hold on;

plot(mdot, T3 - 273.15, 'ob', 'DisplayName', 'T3');
plot(mdot, T4 - 273.15, 'sr', 'DisplayName', 'T4');
plot(mdot, Te_avg- 273.15, 'rs', 'DisplayName', 'Te','MarkerFaceColor', 'r');
plot(mdot, Tc_avg- 273.15, 'ob', 'DisplayName', 'Tc', 'MarkerFaceColor', 'b');
plot(mdot, [Tamb, Tamb, Tamb, Tamb]-273.15, 'k--', 'DisplayName', 'Tamb');
legend();
grid();
title('Figure 1a, Justin Francis');
xlabel('Mass Flow Rate, $\dot m$[kg/s]','Interpreter','latex');
ylabel('Temperature, T[C]');
saveas(gcf, 'Fig1a.png');
%% 1b
figure(2);
hold on;
plot(mdot, qL, 'bo');
plot(mdot, qH, 'rs');
plot(mdot, qloss, 'gx');
plot(mdot, win, 'kd');
legend({'q_L', 'q_H', 'q_{loss}', 'w_{in}'});
grid();
title('Figure 1b, Justin Francis');
xlabel('Mass Flow Rate, $\dot m$[kg/s]','Interpreter','latex');
ylabel('Specific Energy, [kJ/kg]');
saveas(gcf, 'Fig1b.png');

%% 1c
figure(3);
hold on;
plot(mdot,COPR, 'bo');
grid();
title('Figure 1c, Justin Francis');
xlabel('Mass Flow Rate, $\dot m$[kg/s]','Interpreter','latex');
ylabel('COP_R');
saveas(gcf, 'Fig1c.png');

%% 1d
figure(4);
hold on;
yyaxis left;
ylabel('eta_C, [%]');
plot(mdot,etaC, 'bo');
yyaxis right;
plot(mdot,totalPower, 'rs');
grid();
title('Figure 1d, Justin Francis');
xlabel('Mass Flow Rate, $\dot m$[kg/s]','Interpreter','latex');
ylabel('Electrical Power, Win[W]');
legend({'COPR', 'Electrical Power'},'Location','southeast');
saveas(gcf, 'Fig1d.png');
%% 1e
openfig('Ph_Diagram_R134a.fig');
hold on;
maxFlowIdx = 1;
H = [H1(maxFlowIdx), H2(maxFlowIdx), H3(maxFlowIdx), H4(maxFlowIdx), H5(maxFlowIdx)].*10^-3;
P = [P1(maxFlowIdx), P2(maxFlowIdx), P3(maxFlowIdx), P4(maxFlowIdx), P5(maxFlowIdx)].*10^-3;

plot(H, P, 'ro-', 'MarkerFaceColor', 'r');
% plot(H1(maxFlowIdx)*10^-3, P1(maxFlowIdx)*10^-3, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', '1');
% plot(H2(maxFlowIdx)*10^-3, P2(maxFlowIdx)*10^-3, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', '2');
% plot(H3(maxFlowIdx)*10^-3, P3(maxFlowIdx)*10^-3, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', '3');
% plot(H4(maxFlowIdx)*10^-3, P4(maxFlowIdx)*10^-3, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', '4');
% plot(H5(maxFlowIdx)*10^-3, P5(maxFlowIdx)*10^-3, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', '5');
text(H(1), P(1), '1');
text(H(2), P(2), '2');
text(H(3), P(3), '3');
text(H(4), P(4), '4');
text(H(5), P(5), '5');

saveas(gcf, 'Fig1e.png');




