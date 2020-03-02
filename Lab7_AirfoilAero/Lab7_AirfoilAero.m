% Justin Francis
% MEEN 4650, TFES Lab
% Airfoil Aerodynamics

clc; clear; close all;

%% gen vars
n = 17; %number of data files
g = 9.8; %m/s
rho = 1.054; %kg/m3
Patm = 663.4; %mmHg
Tatm = 19; %C
Poutput = 398.144; % Pa
Uinf = ((Poutput*2)/rho)^(0.5); %m/s
c = 4 * 0.0254; %m
s = 12 * 0.0254; %m
Ap = s*c; %m2
x = linspace(1,n,n);
vu = 1.5111e-5; %m2/s
inH20_to_Pa = 248.84; 
Pinf = Poutput;

liftMean = zeros(n,1); 
liftStd = zeros(n,1); 
dragMean = zeros(n,1); 
dragStd = zeros(n,1); 

%% a, find mean and std of lift/drag
for i= [1:n]
    % filename of ith file
    FileName=['./Data/angle_',num2str((i-1)),'.txt'];
    % read in data from the file
    data = importdata(FileName);
    % calculate mean and std
    liftMean(i) = mean(data(:,2)) * g; %N
    liftStd(i) = std(data(:,2)) * g; %N
    dragMean(i) = mean(data(:,3)) * g; %N
    dragStd(i) = std(data(:,3)) * g; %N
end

%% b calc uncertainty
C_L = liftMean/(0.5*rho*Uinf^2*Ap);
C_Lerr = (liftStd/(0.5*rho*Uinf^2*Ap));

C_D = dragMean/(0.5*rho*Uinf^2*Ap);
C_Derr = (dragStd/(0.5*rho*Uinf^2*Ap));



%%%%%%%%%%%%%%%%%%%%% Figs
%% fig 1a
Rec = (Uinf*c)/vu;

figure();
openfig('NACA0012_CL.fig');
hold on;
errorbar(x,C_L.', C_Lerr.', '--bo', 'DisplayName', '1.84e+05');
ylim([0 1.4]);
title('Figure 1a, Justin Francis');
grid();
legend();
saveas(gcf, 'Fig1a.png');

%% fig 1b
figure();
openfig('NACA0012_CD.fig');
hold on;
errorbar(x,abs(C_D.'), C_Derr.', '--bo', 'DisplayName', '1.84e+05');
ylim([0, 0.3]);
title('Figure 1a, Justin Francis');
grid();
legend();
saveas(gcf, 'Fig1b.png');

%% fig 1c
Pout_top = [-2.36, -1.68, -1.27 ,-1.037, -.7, -.55, -.422, -.275, -.129] * inH20_to_Pa;
Pout_bot = [.72 .13 -.14 -.218 -.241 -.247 -.236 -.202 -.154] * inH20_to_Pa;

Pout_top12 = [-1.2 -1.15 -1.1 -1.08 -1.1 -1.05 -1.12 -1.04 -.85] * inH20_to_Pa;
Pout_bot12 = [1.01 .4 .04 -.13 -.23 -.32 -.333 -.45 -.5] * inH20_to_Pa;


Px_top = Pout_top + Pinf;
Px_bot = Pout_bot + Pinf;
Px_top12 = Pout_top12 + Pinf;
Px_bot12 = Pout_bot12 + Pinf;


Cp5top = Pout_top/(.5*rho*Uinf^2);
Cp5bot = Pout_bot/(.5*rho*Uinf^2);
Cp12top = Pout_top12/(0.5*rho*Uinf^2);
Cp12bot = Pout_bot12/(0.5*rho*Uinf^2);

xc = [0.05 .1 .2 .3 .4 .5 .6 .7 .8];

figure();
subplot(1,2,1);
hold on;
plot(xc, Cp5top, 'bo', 'DisplayName', 'Top');
plot(xc, Cp5bot, 'rs', 'DisplayName', 'Bot');
title('Cp for \alpha = 5 [deg]');
xlabel('x/c');
ylabel('C_p');
grid();
legend();

subplot(1,2,2);
hold on;
plot(xc,Cp12top, 'bo', 'DisplayName', 'Top');
plot(xc,Cp12bot, 'rs', 'DisplayName', 'Bot');
title('Cp for \alpha = 12 [deg]');
xlabel('x/c');
ylabel('C_p');
grid();
legend();

saveas(gcf, 'Fig1c.png');







